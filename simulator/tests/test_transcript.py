# python -m unittest tests.test_gene 
import unittest
from simulator.transcript import Transcript

# class GeneSetup(unittest.TestCase):
#     def setUp(self):
#         self.simple_gene = Gene("simulator/tests/testing_data/genes/basic_gene.json")
    
#     def tearDown(self):
#         self.simple_gene = None

class TranscriptCreationTests(unittest.TestCase): 
    def test_normal_gene(self):
        transcript = Transcript("GALK2", 0)
        print transcript

    def test_assert_raises(self):
        self.assertRaises(Transcript("random", 0))

class TrascriptMethodTests(unittest.TestCase):
    def test_get_seq(self):
        transcript = Transcript("SOX9", 0)
        seq = transcript.get_seq()
        self.assertEquals(seq, "GGAGAGCCGAAAGCGGAGCTCGAAACTGACTGGAAACTTCAGTGGCGCGGAGACTCGCCAGTTTCAACCCCGGAAACTTTTCTTTGCAGGAGGAGAAGAGAAGGGGTGCAAGCGCCCCCACTTTTGCTCTTTTTCCTCCCCTCCTCCTCCTCTCCAATTCGCCTCCCCCCACTTGGAGCGGGCAGCTGTGAACTGGCCACCCCGCGCCTTCCTAAGTGCTCGCCGCGGTAGCCGGCCGACGCGCCAGCTTCCCCGGGAGCCGCTTGCTCCGCATCCGGGCAGCCGAGGGGAGAGGAGCCCGCGCCTCGAGTCCCCGAGCCGCCGCGGCTTCTCGCCTTTCCCGGCCACCAGCCCCCTGCCCCGGGCCCGCGTATGAATCTCCTGGACCCCTTCATGAAGATGACCGACGAGCAGGAGAAGGGCCTGTCCGGCGCCCCCAGCCCCACCATGTCCGAGGACTCCGCGGGCTCGCCCTGCCCGTCGGGCTCCGGCTCGGACACCGAGAACACGCGGCCCCAGGAGAACACGTTCCCCAAGGGCGAGCCCGATCTGAAGAAGGAGAGCGAGGAGGACAAGTTCCCCGTGTGCATCCGCGAGGCGGTCAGCCAGGTGCTCAAAGGCTACGACTGGACGCTGGTGCCCATGCCGGTGCGCGTCAACGGCTCCAGCAAGAACAAGCCGCACGTCAAGCGGCCCATGAACGCCTTCATGGTGTGGGCGCAGGCGGCGCGCAGGAAGCTCGCGGACCAGTACCCGCACTTGCACAACGCCGAGCTCAGCAAGACGCTGGGCAAGCTCTGGAGgtaggacccggcgggggcggcgcggcagggtgggcatcgcggcggctgggggcgctggtcagggctgatttgccccgccccgcctcccatcgcccgggagttgccgttccgggagccggcgggatggggttgggagtgggaatggggtgtaactgtggctcagagtttgacaaagttcttgggctgctcgcggggacgcggaggaggggggtggtaagtggaagaggtgagggaggtagctggaggatggacgaagactggtgggagacggaaggagggggctgccagcctgctctccagtcgcctggaagctcaatcggggcggggaagtgaaacttgcctccctcctacccggcctcttaaaactgcactctctcgtgcagccccactgtccacggagatggggcaagggagaaaccgaggttggaggagacccttggcaggaactgggaggcgggaggagggaggctactggaaataggtgggagtgtatggtggggggtgagaattggggaccttcttgcagcttaagtaatttgggggaaagttttcaaagggggttggggttgggggcggtaagtcgagcagcaaaggcgtttagggggcagcaccgggagtcgttttcatctccagcgtttccaaaatagaaatagaaggggaggggagggagggggcggggagtgaccgctcaggtcagactgcaataacttatttatttatttatttttaagaaaagttatgagctgtggttgcaggcaggagggaagatggagttgtgtgcagaggaagccgagtggtctgggtcgccgcctcctccccgccgacctgacagtttggcggatttcactgacccctctccctctttttctctgtgccccccgccccgccccgagcagACTTCTGAACGAGAGCGAGAAGCGGCCCTTCGTGGAGGAGGCGGAGCGGCTGCGCGTGCAGCACAAGAAGGACCACCCGGATTACAAGTACCAGCCGCGGCGGAGGAAGTCGGTGAAGAACGGGCAGGCGGAGGCAGAGGAGGCCACGGAGCAGACGCACATCTCCCCCAACGCCATCTTCAAGGCGCTGCAGGCCGACTCGCCACACTCCTCCTCCGGCATGAGCGAGGTGCACTCCCCCGGCGAGCACTCGGgtgagtcgcccctcgaccccaccggacaagctatctccgtcccgcctggcacaccccctgccctccgcctgggagattcttcgtggggactttatgcttcccgggagggacacactgccctttgcgcccgtcccgctcccctctctacccagagcctaagaggcatccaaacaacacacacacaaacacacacaccccaactcaatcccagcatccgaagagattaacttttttattgggaggtaaaatgcccttaacagccttacaagacctctcccttcttctctgctcccccaccccaaaagcacacacagggctcttacacaagtagcaattaggtcttccggaccctccgggccccagaccctcccctgataaaagggggctgtccagtgtgtaccggcgggttaatcattgggcgacttatctccggtgcagcgcgcctcttgcgcgggtgcgggcccttattacactttagcagcgagggagggtccccggagggtgcctaagactagggcgtctgcacagcccttgttgattttctcgtgcttgttcttttattgtccacagGGCAATCCCAGGGCCCACCGACCCCACCCACCACCCCCAAAACCGACGTGCAGCCGGGCAAGGCTGACCTGAAGCGAGAGGGGCGCCCCTTGCCAGAGGGGGGCAGACAGCCCCCTATCGACTTCCGCGACGTGGACATCGGCGAGCTGAGCAGCGACGTCATCTCCAACATCGAGACCTTCGATGTCAACGAGTTTGACCAGTACCTGCCGCCCAACGGCCACCCGGGGGTGCCGGCCACGCACGGCCAGGTCACCTACACGGGCAGCTACGGCATCAGCAGCACCGCGGCCACCCCGGCGAGCGCGGGCCACGTGTGGATGTCCAAGCAGCAGGCGCCGCCGCCACCCCCGCAGCAGCCCCCACAGGCCCCGCCGGCCCCGCAGGCGCCCCCGCAGCCGCAGGCGGCGCCCCCACAGCAGCCGGCGGCACCCCCGCAGCAGCCACAGGCGCACACGCTGACCACGCTGAGCAGCGAGCCGGGCCAGTCCCAGCGAACGCACATCAAGACGGAGCAGCTGAGCCCCAGCCACTACAGCGAGCAGCAGCAGCACTCGCCCCAACAGATCGCCTACAGCCCCTTCAACCTCCCACACTACAGCCCCTCCTACCCGCCCATCACCCGCTCACAGTACGACTACACCGACCACCAGAACTCCAGCTCCTACTACAGCCACGCGGCAGGCCAGGGCACCGGCCTCTACTCCACCTTCACCTACATGAACCCCGCTCAGCGCCCCATGTACACCCCCATCGCCGACACCTCTGGGGTCCCTTCCATCCCGCAGACCCACAGCCCCCAGCACTGGGAACAACCCGTCTACACACAGCTCACTCGACCTTGAGGAGGCCTCCCACGAAGGGCGAAGATGGCCGAGATGATCCTAAAAATAACCGAAGAAAGAGAGGACCAACCAGAATTCCCTTTGGACATTTGTGTTTTTTTGTTTTTTTATTTTGTTTTGTTTTTTCTTCTTCTTCTTCTTCCTTAAAGACATTTAAGCTAAAGGCAACTCGTACCCAAATTTCCAAGACACAAACATGACCTATCCAAGCGCATTACCCACTTGTGGCCAATCAGTGGCCAGGCCAACCTTGGCTAAATGGAGCAGCGAAATCAACGAGAAACTGGACTTTTTAAACCCTCTTCAGAGCAAGCGTGGAGGATGATGGAGAATCGTGTGATCAGTGTGCTAAATCTCTCTGCCTGTTTGGACTTTGTAATTATTTTTTTAGCAGTAATTAAAGAAAAAAGTCCTCTGTGAGGAATATTCTCTATTTTAAATATTTTTAGTATGTACTGTGTATGATTCATTACCATTTTGAGGGGATTTATACATATTTTTAGATAAAATTAAATGCTCTTATTTTTCCAACAGCTAAACTACTCTTAGTTGAACAGTGTGCCCTAGCTTTTCTTGCAACCAGAGTATTTTTGTACAGATTTGCTTTCTCTTACAAAAAGAAAAAAAAAATCCTGTTGTATTAACATTTAAAAACAGAATTGTGTTATGTGATCAGTTTTGGGGGTTAACTTTGCTTAATTCCTCAGGCTTTGCGATTTAAGGAGGAGCTGCCTTAAAAAAAAATAAAGGCCTTATTTTGCAATTATGGGAGTAAACAATAGTCTAGAGAAGCATTTGGTAAGCTTTATCATATATATATTTTTTAAAGAAGAGAAAAACACCTTGAGCCTTAAAACGGTGCTGCTGGGAAACATTTGCACTCTTTTAGTGCATTTCCTCCTGCCTTTGCTTGTTCACTGCAGTCTTAAGAAAGAGGTAAAAGGCAAGCAAAGGAGATGAAATCTGTTCTGGGAATGTTTCAGCAGCCAATAAGTGCCCGAGCACACTGCCCCCGGTTGCCTGCCTGGGCCCCATGTGGAAGGCAGATGCCTGCTCGCTCTGTCACCTGTGCCTCTCAGAACACCAGCAGTTAACCTTCAAGACATTCCACTTGCTAAAATTATTTATTTTGTAAGGAGAGGTTTTAATTAAAACAAAAAAAAATTCTTTTTTTTTTTTTTTTCCAATTTTACCTTCTTTAAAATAGGTTGTTGGAGCTTTCCTCAAAGGGTATGGTCATCTGTTGTTAAATTATGTTCTTAACTGTAACCAGTTTTTTTTTATTTATCTCTTTAATCTTTTTTTATTATTAAAAGCAAGTTTCTTTGTATTCCTCACCCTAGATTTGTATAAATGCCTTTTTGTCCATCCCTTTTTTCTTTGTTGTTTTTGTTGAAAACAAACTGGAAACTTGTTTCTTTTTTTGTATAAATGAGAGATTGCAAATGTAGTGTATCACTGAGTCATTTGCAGTGTTTTCTGCCACAGACCTTTGGGCTGCCTTATATTGTGTGTGTGTGTGGGTGTGTGTGTGTTTTGACACAAAAACAATGCAAGCATGTGTCATCCATATTTCTCTGCATCTTCTCTTGGAGTGAGGGAGGCTACCTGGAGGGGATCAGCCCACTGACAGACCTTAATCTTAATTACTGCTGTGGCTAGAGAGTTTGAGGATTGCTTTTTAAAAAAGACAGCAAACTTTTTTTTTTATTTAAAAAAAGATATATTAACAGTTTTAGAAGTCAGTAGAATAAAATCTTAAAGCACTCATAATATGGCATCCTTCAATTTCTGTATAAAAGCAGATCTTTTTAAAAAGATACTTCTGTAACTTAAGAAACCTGGCATTTAAATCATATTTTGTCTTTAGGTAAAAGCTTTGGTTTGTGTTCGTGTTTTGTTTGTTTCACTTGTTTCCCTCCCAGCCCCAAACCTTTTGTTCTCTCCGTGAAACTTACCTTTCCCTTTTTCTTTCTCTTTTTTTTTTTTGTATATTATTGTTTACAATAAATATACATTGCATTAAAAAGAA".upper())

    def test_get_exons(self):
        transcript = Transcript("SOX9", 0)
        exons = transcript.get_exons()
        exon_seq = transcript.get_seq_from_pos(exons)
        self.assertEquals(exon_seq, "GGAGAGCCGAAAGCGGAGCTCGAAACTGACTGGAAACTTCAGTGGCGCGGAGACTCGCCAGTTTCAACCCCGGAAACTTTTCTTTGCAGGAGGAGAAGAGAAGGGGTGCAAGCGCCCCCACTTTTGCTCTTTTTCCTCCCCTCCTCCTCCTCTCCAATTCGCCTCCCCCCACTTGGAGCGGGCAGCTGTGAACTGGCCACCCCGCGCCTTCCTAAGTGCTCGCCGCGGTAGCCGGCCGACGCGCCAGCTTCCCCGGGAGCCGCTTGCTCCGCATCCGGGCAGCCGAGGGGAGAGGAGCCCGCGCCTCGAGTCCCCGAGCCGCCGCGGCTTCTCGCCTTTCCCGGCCACCAGCCCCCTGCCCCGGGCCCGCGTATGAATCTCCTGGACCCCTTCATGAAGATGACCGACGAGCAGGAGAAGGGCCTGTCCGGCGCCCCCAGCCCCACCATGTCCGAGGACTCCGCGGGCTCGCCCTGCCCGTCGGGCTCCGGCTCGGACACCGAGAACACGCGGCCCCAGGAGAACACGTTCCCCAAGGGCGAGCCCGATCTGAAGAAGGAGAGCGAGGAGGACAAGTTCCCCGTGTGCATCCGCGAGGCGGTCAGCCAGGTGCTCAAAGGCTACGACTGGACGCTGGTGCCCATGCCGGTGCGCGTCAACGGCTCCAGCAAGAACAAGCCGCACGTCAAGCGGCCCATGAACGCCTTCATGGTGTGGGCGCAGGCGGCGCGCAGGAAGCTCGCGGACCAGTACCCGCACTTGCACAACGCCGAGCTCAGCAAGACGCTGGGCAAGCTCTGGAGACTTCTGAACGAGAGCGAGAAGCGGCCCTTCGTGGAGGAGGCGGAGCGGCTGCGCGTGCAGCACAAGAAGGACCACCCGGATTACAAGTACCAGCCGCGGCGGAGGAAGTCGGTGAAGAACGGGCAGGCGGAGGCAGAGGAGGCCACGGAGCAGACGCACATCTCCCCCAACGCCATCTTCAAGGCGCTGCAGGCCGACTCGCCACACTCCTCCTCCGGCATGAGCGAGGTGCACTCCCCCGGCGAGCACTCGGGGCAATCCCAGGGCCCACCGACCCCACCCACCACCCCCAAAACCGACGTGCAGCCGGGCAAGGCTGACCTGAAGCGAGAGGGGCGCCCCTTGCCAGAGGGGGGCAGACAGCCCCCTATCGACTTCCGCGACGTGGACATCGGCGAGCTGAGCAGCGACGTCATCTCCAACATCGAGACCTTCGATGTCAACGAGTTTGACCAGTACCTGCCGCCCAACGGCCACCCGGGGGTGCCGGCCACGCACGGCCAGGTCACCTACACGGGCAGCTACGGCATCAGCAGCACCGCGGCCACCCCGGCGAGCGCGGGCCACGTGTGGATGTCCAAGCAGCAGGCGCCGCCGCCACCCCCGCAGCAGCCCCCACAGGCCCCGCCGGCCCCGCAGGCGCCCCCGCAGCCGCAGGCGGCGCCCCCACAGCAGCCGGCGGCACCCCCGCAGCAGCCACAGGCGCACACGCTGACCACGCTGAGCAGCGAGCCGGGCCAGTCCCAGCGAACGCACATCAAGACGGAGCAGCTGAGCCCCAGCCACTACAGCGAGCAGCAGCAGCACTCGCCCCAACAGATCGCCTACAGCCCCTTCAACCTCCCACACTACAGCCCCTCCTACCCGCCCATCACCCGCTCACAGTACGACTACACCGACCACCAGAACTCCAGCTCCTACTACAGCCACGCGGCAGGCCAGGGCACCGGCCTCTACTCCACCTTCACCTACATGAACCCCGCTCAGCGCCCCATGTACACCCCCATCGCCGACACCTCTGGGGTCCCTTCCATCCCGCAGACCCACAGCCCCCAGCACTGGGAACAACCCGTCTACACACAGCTCACTCGACCTTGAGGAGGCCTCCCACGAAGGGCGAAGATGGCCGAGATGATCCTAAAAATAACCGAAGAAAGAGAGGACCAACCAGAATTCCCTTTGGACATTTGTGTTTTTTTGTTTTTTTATTTTGTTTTGTTTTTTCTTCTTCTTCTTCTTCCTTAAAGACATTTAAGCTAAAGGCAACTCGTACCCAAATTTCCAAGACACAAACATGACCTATCCAAGCGCATTACCCACTTGTGGCCAATCAGTGGCCAGGCCAACCTTGGCTAAATGGAGCAGCGAAATCAACGAGAAACTGGACTTTTTAAACCCTCTTCAGAGCAAGCGTGGAGGATGATGGAGAATCGTGTGATCAGTGTGCTAAATCTCTCTGCCTGTTTGGACTTTGTAATTATTTTTTTAGCAGTAATTAAAGAAAAAAGTCCTCTGTGAGGAATATTCTCTATTTTAAATATTTTTAGTATGTACTGTGTATGATTCATTACCATTTTGAGGGGATTTATACATATTTTTAGATAAAATTAAATGCTCTTATTTTTCCAACAGCTAAACTACTCTTAGTTGAACAGTGTGCCCTAGCTTTTCTTGCAACCAGAGTATTTTTGTACAGATTTGCTTTCTCTTACAAAAAGAAAAAAAAAATCCTGTTGTATTAACATTTAAAAACAGAATTGTGTTATGTGATCAGTTTTGGGGGTTAACTTTGCTTAATTCCTCAGGCTTTGCGATTTAAGGAGGAGCTGCCTTAAAAAAAAATAAAGGCCTTATTTTGCAATTATGGGAGTAAACAATAGTCTAGAGAAGCATTTGGTAAGCTTTATCATATATATATTTTTTAAAGAAGAGAAAAACACCTTGAGCCTTAAAACGGTGCTGCTGGGAAACATTTGCACTCTTTTAGTGCATTTCCTCCTGCCTTTGCTTGTTCACTGCAGTCTTAAGAAAGAGGTAAAAGGCAAGCAAAGGAGATGAAATCTGTTCTGGGAATGTTTCAGCAGCCAATAAGTGCCCGAGCACACTGCCCCCGGTTGCCTGCCTGGGCCCCATGTGGAAGGCAGATGCCTGCTCGCTCTGTCACCTGTGCCTCTCAGAACACCAGCAGTTAACCTTCAAGACATTCCACTTGCTAAAATTATTTATTTTGTAAGGAGAGGTTTTAATTAAAACAAAAAAAAATTCTTTTTTTTTTTTTTTTCCAATTTTACCTTCTTTAAAATAGGTTGTTGGAGCTTTCCTCAAAGGGTATGGTCATCTGTTGTTAAATTATGTTCTTAACTGTAACCAGTTTTTTTTTATTTATCTCTTTAATCTTTTTTTATTATTAAAAGCAAGTTTCTTTGTATTCCTCACCCTAGATTTGTATAAATGCCTTTTTGTCCATCCCTTTTTTCTTTGTTGTTTTTGTTGAAAACAAACTGGAAACTTGTTTCTTTTTTTGTATAAATGAGAGATTGCAAATGTAGTGTATCACTGAGTCATTTGCAGTGTTTTCTGCCACAGACCTTTGGGCTGCCTTATATTGTGTGTGTGTGTGGGTGTGTGTGTGTTTTGACACAAAAACAATGCAAGCATGTGTCATCCATATTTCTCTGCATCTTCTCTTGGAGTGAGGGAGGCTACCTGGAGGGGATCAGCCCACTGACAGACCTTAATCTTAATTACTGCTGTGGCTAGAGAGTTTGAGGATTGCTTTTTAAAAAAGACAGCAAACTTTTTTTTTTATTTAAAAAAAGATATATTAACAGTTTTAGAAGTCAGTAGAATAAAATCTTAAAGCACTCATAATATGGCATCCTTCAATTTCTGTATAAAAGCAGATCTTTTTAAAAAGATACTTCTGTAACTTAAGAAACCTGGCATTTAAATCATATTTTGTCTTTAGGTAAAAGCTTTGGTTTGTGTTCGTGTTTTGTTTGTTTCACTTGTTTCCCTCCCAGCCCCAAACCTTTTGTTCTCTCCGTGAAACTTACCTTTCCCTTTTTCTTTCTCTTTTTTTTTTTTGTATATTATTGTTTACAATAAATATACATTGCATTAAAAAGAA")

    def test_get_exons_negative(self):
        transcript = Transcript("TOR1A", 0)
        exons = transcript.get_exons()
        exon_seq = transcript.get_seq_from_pos(exons)
        self.assertEquals(exon_seq, "GAACCGGAAGCGTGGGTCTGGCGGCTGCACCGGTTCGCGGTCGGCGCGAGAACAAGCAGGGTGGCGCGGGTCCGGGCATGAAGCTGGGCCGGGCCGTGCTGGGCCTGCTGCTGCTGGCGCCGTCCGTGGTGCAGGCGGTGGAGCCCATCAGCCTGGGACTGGCCCTGGCCGGCGTCCTCACCGGCTACATCTACCCGCGTCTCTACTGCCTCTTCGCCGAGTGCTGCGGGCAGAAGCGGAGCCTTAGCCGGGAGGCACTGCAGAAGGATCTGGACGACAACCTCTTTGGACAGCATCTTGCAAAGAAAATCATCTTAAATGCCGTGTTTGGTTTCATAAACAACCCAAAGCCCAAGAAACCTCTCACGCTCTCCCTGCACGGGTGGACAGGCACCGGCAAAAATTTCGTCAGCAAGATCATCGCAGAGAATATTTACGAGGGTGGTCTGAACAGTGACTATGTCCACCTGTTTGTGGCCACATTGCACTTTCCACATGCTTCAAACATCACCTTGTACAAGGATCAGTTACAGTTGTGGATTCGAGGCAACGTGAGTGCCTGTGCGAGGTCCATCTTCATATTTGATGAAATGGATAAGATGCATGCAGGCCTCATAGATGCCATCAAGCCTTTCCTCGACTATTATGACCTGGTGGATGGGGTCTCCTACCAGAAAGCCATGTTCATATTTCTCAGCAATGCTGGAGCAGAAAGGATCACAGATGTGGCTTTGGATTTCTGGAGGAGTGGAAAGCAGAGGGAAGACATCAAGCTCAAAGACATTGAACACGCGTTGTCTGTGTCGGTTTTCAATAACAAGAACAGTGGCTTCTGGCACAGCAGCTTAATTGACCGGAACCTCATTGATTATTTTGTTCCCTTCCTCCCCCTGGAATACAAACACCTAAAAATGTGTATCCGAGTGGAAATGCAGTCCCGAGGCTATGAAATTGATGAAGACATTGTAAGCAGAGTGGCTGAGGAGATGACATTTTTCCCCAAAGAGGAGAGAGTTTTCTCAGATAAAGGCTGCAAAACGGTGTTCACCAAGTTAGATTATTACTACGATGATTGACAGTCATGATTGGCAGCCGGAGTCACTGCCTGGAGTTGGAAAAGAAACAACACTCAGTCCTTCCACACTTCCACCCCCAGCTCCTTTCCCTGGAAGAGGAATCCAGTGAATGTTCCTGTTTGATGTGACAGGAATTCTCCCTGGCATTGTTTCCACCCCCTGGTGCCTGCAGGCCACCCAGGGACCACGGGCGAGGACGTGAAGCCTCCCGAACACGCACAGAAGGAAGGAGCCAGCTCCCAGCCCACTCATCGCAGGGCTCATGATTTTTTACAAATTATGTTTTAATTCCAAGTGTTTCTGTTTCAAGGAAGGATGAATAAGTTTTATTGAAAATGTGGTAACTTTATTTAAAATGATTTTTAACATTATGAGAGACTGCTCAGATTCTAAGTTGTTGGCCTTGTGTGTGTGTTTTTTTTTAAGTTCTCATCATTATTACATAGACTGTGATGTATCTTTACTGGAAATGAGCCCAAGCACACATGCATGGCATTTGTTCCACAGGAGGGCATCCCTGGGGATGTGGCTGGAGCATGAGCCAGCTCTGTCCCAGGATGGTCCCAGCGGATGCTGCCAGGGGCAGTGAAGTGTTTAGGTGAAGGACAAGTAGGTAAGAGGACGCCTTCAGGCACCACAGATAAGCCTGAAACAGCCTCTCCAAGGGTTTTCACCTTAGCAACAATGGGAGCTGTGGGAGTGATTTTGGCCACACTGTCAACATTTGTTAGAACCAGTCTTTTGAAAGAAAAGTATTTCCAACTTGTCACTTGCCAGTCACTCCGTTTTGCAAAAGGTGGCCCTTCACTGTCCATTCCAAATAGCCCACACGTGCTCTCTGCTGGATTCTAAATTATGTGAATTTTGCCATATTAAATCTTCCTCATTTATACTATTATTTGTTACGTTCAATCAGAATCCCCGAAACCTCCTATAAAGCTTAGCTGCCCCTTCTGAGGATGCTGAGAACGGTGTCTTTCTTTATAAATGCAAATGGCTACCGTTTTACAATAAAATTTTGCATGTGCCA")

    def test_get_coding_exons(self):
        transcript = Transcript("SOX9", 0)
        exons = transcript.get_coding()
        exon_seq = transcript.get_seq_from_pos(exons)
        self.assertEquals(exon_seq, "ATGAATCTCCTGGACCCCTTCATGAAGATGACCGACGAGCAGGAGAAGGGCCTGTCCGGCGCCCCCAGCCCCACCATGTCCGAGGACTCCGCGGGCTCGCCCTGCCCGTCGGGCTCCGGCTCGGACACCGAGAACACGCGGCCCCAGGAGAACACGTTCCCCAAGGGCGAGCCCGATCTGAAGAAGGAGAGCGAGGAGGACAAGTTCCCCGTGTGCATCCGCGAGGCGGTCAGCCAGGTGCTCAAAGGCTACGACTGGACGCTGGTGCCCATGCCGGTGCGCGTCAACGGCTCCAGCAAGAACAAGCCGCACGTCAAGCGGCCCATGAACGCCTTCATGGTGTGGGCGCAGGCGGCGCGCAGGAAGCTCGCGGACCAGTACCCGCACTTGCACAACGCCGAGCTCAGCAAGACGCTGGGCAAGCTCTGGAGACTTCTGAACGAGAGCGAGAAGCGGCCCTTCGTGGAGGAGGCGGAGCGGCTGCGCGTGCAGCACAAGAAGGACCACCCGGATTACAAGTACCAGCCGCGGCGGAGGAAGTCGGTGAAGAACGGGCAGGCGGAGGCAGAGGAGGCCACGGAGCAGACGCACATCTCCCCCAACGCCATCTTCAAGGCGCTGCAGGCCGACTCGCCACACTCCTCCTCCGGCATGAGCGAGGTGCACTCCCCCGGCGAGCACTCGGGGCAATCCCAGGGCCCACCGACCCCACCCACCACCCCCAAAACCGACGTGCAGCCGGGCAAGGCTGACCTGAAGCGAGAGGGGCGCCCCTTGCCAGAGGGGGGCAGACAGCCCCCTATCGACTTCCGCGACGTGGACATCGGCGAGCTGAGCAGCGACGTCATCTCCAACATCGAGACCTTCGATGTCAACGAGTTTGACCAGTACCTGCCGCCCAACGGCCACCCGGGGGTGCCGGCCACGCACGGCCAGGTCACCTACACGGGCAGCTACGGCATCAGCAGCACCGCGGCCACCCCGGCGAGCGCGGGCCACGTGTGGATGTCCAAGCAGCAGGCGCCGCCGCCACCCCCGCAGCAGCCCCCACAGGCCCCGCCGGCCCCGCAGGCGCCCCCGCAGCCGCAGGCGGCGCCCCCACAGCAGCCGGCGGCACCCCCGCAGCAGCCACAGGCGCACACGCTGACCACGCTGAGCAGCGAGCCGGGCCAGTCCCAGCGAACGCACATCAAGACGGAGCAGCTGAGCCCCAGCCACTACAGCGAGCAGCAGCAGCACTCGCCCCAACAGATCGCCTACAGCCCCTTCAACCTCCCACACTACAGCCCCTCCTACCCGCCCATCACCCGCTCACAGTACGACTACACCGACCACCAGAACTCCAGCTCCTACTACAGCCACGCGGCAGGCCAGGGCACCGGCCTCTACTCCACCTTCACCTACATGAACCCCGCTCAGCGCCCCATGTACACCCCCATCGCCGACACCTCTGGGGTCCCTTCCATCCCGCAGACCCACAGCCCCCAGCACTGGGAACAACCCGTCTACACACAGCTCACTCGACCTTGA")
    
    def test_noncoding_gene(self):
        transcript = Transcript("PRDM15", 0)
        exons = transcript.get_coding()
        exon_seq = transcript.get_seq_from_pos(exons)
        self.assertEquals(exon_seq, "")


    def test_get_introns(self):
        transcript = Transcript("SOX9", 0)
        introns = transcript.get_introns()
        intron_seq = transcript.get_seq_from_pos(introns)
        self.assertEquals(intron_seq, "gtaggacccggcgggggcggcgcggcagggtgggcatcgcggcggctgggggcgctggtcagggctgatttgccccgccccgcctcccatcgcccgggagttgccgttccgggagccggcgggatggggttgggagtgggaatggggtgtaactgtggctcagagtttgacaaagttcttgggctgctcgcggggacgcggaggaggggggtggtaagtggaagaggtgagggaggtagctggaggatggacgaagactggtgggagacggaaggagggggctgccagcctgctctccagtcgcctggaagctcaatcggggcggggaagtgaaacttgcctccctcctacccggcctcttaaaactgcactctctcgtgcagccccactgtccacggagatggggcaagggagaaaccgaggttggaggagacccttggcaggaactgggaggcgggaggagggaggctactggaaataggtgggagtgtatggtggggggtgagaattggggaccttcttgcagcttaagtaatttgggggaaagttttcaaagggggttggggttgggggcggtaagtcgagcagcaaaggcgtttagggggcagcaccgggagtcgttttcatctccagcgtttccaaaatagaaatagaaggggaggggagggagggggcggggagtgaccgctcaggtcagactgcaataacttatttatttatttatttttaagaaaagttatgagctgtggttgcaggcaggagggaagatggagttgtgtgcagaggaagccgagtggtctgggtcgccgcctcctccccgccgacctgacagtttggcggatttcactgacccctctccctctttttctctgtgccccccgccccgccccgagcaggtgagtcgcccctcgaccccaccggacaagctatctccgtcccgcctggcacaccccctgccctccgcctgggagattcttcgtggggactttatgcttcccgggagggacacactgccctttgcgcccgtcccgctcccctctctacccagagcctaagaggcatccaaacaacacacacacaaacacacacaccccaactcaatcccagcatccgaagagattaacttttttattgggaggtaaaatgcccttaacagccttacaagacctctcccttcttctctgctcccccaccccaaaagcacacacagggctcttacacaagtagcaattaggtcttccggaccctccgggccccagaccctcccctgataaaagggggctgtccagtgtgtaccggcgggttaatcattgggcgacttatctccggtgcagcgcgcctcttgcgcgggtgcgggcccttattacactttagcagcgagggagggtccccggagggtgcctaagactagggcgtctgcacagcccttgttgattttctcgtgcttgttcttttattgtccacag".upper())

    def test_get_utrs_positive_strand(self):
        transcript = Transcript("SOX9", 0)
        both_utr = transcript.get_utr("both")
        both_utr = transcript.get_seq_from_pos(both_utr)
        self.assertEquals(both_utr, "GGAGAGCCGAAAGCGGAGCTCGAAACTGACTGGAAACTTCAGTGGCGCGGAGACTCGCCAGTTTCAACCCCGGAAACTTTTCTTTGCAGGAGGAGAAGAGAAGGGGTGCAAGCGCCCCCACTTTTGCTCTTTTTCCTCCCCTCCTCCTCCTCTCCAATTCGCCTCCCCCCACTTGGAGCGGGCAGCTGTGAACTGGCCACCCCGCGCCTTCCTAAGTGCTCGCCGCGGTAGCCGGCCGACGCGCCAGCTTCCCCGGGAGCCGCTTGCTCCGCATCCGGGCAGCCGAGGGGAGAGGAGCCCGCGCCTCGAGTCCCCGAGCCGCCGCGGCTTCTCGCCTTTCCCGGCCACCAGCCCCCTGCCCCGGGCCCGCGTGGAGGCCTCCCACGAAGGGCGAAGATGGCCGAGATGATCCTAAAAATAACCGAAGAAAGAGAGGACCAACCAGAATTCCCTTTGGACATTTGTGTTTTTTTGTTTTTTTATTTTGTTTTGTTTTTTCTTCTTCTTCTTCTTCCTTAAAGACATTTAAGCTAAAGGCAACTCGTACCCAAATTTCCAAGACACAAACATGACCTATCCAAGCGCATTACCCACTTGTGGCCAATCAGTGGCCAGGCCAACCTTGGCTAAATGGAGCAGCGAAATCAACGAGAAACTGGACTTTTTAAACCCTCTTCAGAGCAAGCGTGGAGGATGATGGAGAATCGTGTGATCAGTGTGCTAAATCTCTCTGCCTGTTTGGACTTTGTAATTATTTTTTTAGCAGTAATTAAAGAAAAAAGTCCTCTGTGAGGAATATTCTCTATTTTAAATATTTTTAGTATGTACTGTGTATGATTCATTACCATTTTGAGGGGATTTATACATATTTTTAGATAAAATTAAATGCTCTTATTTTTCCAACAGCTAAACTACTCTTAGTTGAACAGTGTGCCCTAGCTTTTCTTGCAACCAGAGTATTTTTGTACAGATTTGCTTTCTCTTACAAAAAGAAAAAAAAAATCCTGTTGTATTAACATTTAAAAACAGAATTGTGTTATGTGATCAGTTTTGGGGGTTAACTTTGCTTAATTCCTCAGGCTTTGCGATTTAAGGAGGAGCTGCCTTAAAAAAAAATAAAGGCCTTATTTTGCAATTATGGGAGTAAACAATAGTCTAGAGAAGCATTTGGTAAGCTTTATCATATATATATTTTTTAAAGAAGAGAAAAACACCTTGAGCCTTAAAACGGTGCTGCTGGGAAACATTTGCACTCTTTTAGTGCATTTCCTCCTGCCTTTGCTTGTTCACTGCAGTCTTAAGAAAGAGGTAAAAGGCAAGCAAAGGAGATGAAATCTGTTCTGGGAATGTTTCAGCAGCCAATAAGTGCCCGAGCACACTGCCCCCGGTTGCCTGCCTGGGCCCCATGTGGAAGGCAGATGCCTGCTCGCTCTGTCACCTGTGCCTCTCAGAACACCAGCAGTTAACCTTCAAGACATTCCACTTGCTAAAATTATTTATTTTGTAAGGAGAGGTTTTAATTAAAACAAAAAAAAATTCTTTTTTTTTTTTTTTTCCAATTTTACCTTCTTTAAAATAGGTTGTTGGAGCTTTCCTCAAAGGGTATGGTCATCTGTTGTTAAATTATGTTCTTAACTGTAACCAGTTTTTTTTTATTTATCTCTTTAATCTTTTTTTATTATTAAAAGCAAGTTTCTTTGTATTCCTCACCCTAGATTTGTATAAATGCCTTTTTGTCCATCCCTTTTTTCTTTGTTGTTTTTGTTGAAAACAAACTGGAAACTTGTTTCTTTTTTTGTATAAATGAGAGATTGCAAATGTAGTGTATCACTGAGTCATTTGCAGTGTTTTCTGCCACAGACCTTTGGGCTGCCTTATATTGTGTGTGTGTGTGGGTGTGTGTGTGTTTTGACACAAAAACAATGCAAGCATGTGTCATCCATATTTCTCTGCATCTTCTCTTGGAGTGAGGGAGGCTACCTGGAGGGGATCAGCCCACTGACAGACCTTAATCTTAATTACTGCTGTGGCTAGAGAGTTTGAGGATTGCTTTTTAAAAAAGACAGCAAACTTTTTTTTTTATTTAAAAAAAGATATATTAACAGTTTTAGAAGTCAGTAGAATAAAATCTTAAAGCACTCATAATATGGCATCCTTCAATTTCTGTATAAAAGCAGATCTTTTTAAAAAGATACTTCTGTAACTTAAGAAACCTGGCATTTAAATCATATTTTGTCTTTAGGTAAAAGCTTTGGTTTGTGTTCGTGTTTTGTTTGTTTCACTTGTTTCCCTCCCAGCCCCAAACCTTTTGTTCTCTCCGTGAAACTTACCTTTCCCTTTTTCTTTCTCTTTTTTTTTTTTGTATATTATTGTTTACAATAAATATACATTGCATTAAAAAGAA") 
        prime_5 = transcript.get_utr("5_prime")
        prime_5 = transcript.get_seq_from_pos(prime_5)
        self.assertEquals(prime_5,"GGAGAGCCGAAAGCGGAGCTCGAAACTGACTGGAAACTTCAGTGGCGCGGAGACTCGCCAGTTTCAACCCCGGAAACTTTTCTTTGCAGGAGGAGAAGAGAAGGGGTGCAAGCGCCCCCACTTTTGCTCTTTTTCCTCCCCTCCTCCTCCTCTCCAATTCGCCTCCCCCCACTTGGAGCGGGCAGCTGTGAACTGGCCACCCCGCGCCTTCCTAAGTGCTCGCCGCGGTAGCCGGCCGACGCGCCAGCTTCCCCGGGAGCCGCTTGCTCCGCATCCGGGCAGCCGAGGGGAGAGGAGCCCGCGCCTCGAGTCCCCGAGCCGCCGCGGCTTCTCGCCTTTCCCGGCCACCAGCCCCCTGCCCCGGGCCCGCGT") 
        prime_3 = transcript.get_utr("3_prime")
        prime_3 = transcript.get_seq_from_pos(prime_3)
        self.assertEquals(prime_3,"GGAGGCCTCCCACGAAGGGCGAAGATGGCCGAGATGATCCTAAAAATAACCGAAGAAAGAGAGGACCAACCAGAATTCCCTTTGGACATTTGTGTTTTTTTGTTTTTTTATTTTGTTTTGTTTTTTCTTCTTCTTCTTCTTCCTTAAAGACATTTAAGCTAAAGGCAACTCGTACCCAAATTTCCAAGACACAAACATGACCTATCCAAGCGCATTACCCACTTGTGGCCAATCAGTGGCCAGGCCAACCTTGGCTAAATGGAGCAGCGAAATCAACGAGAAACTGGACTTTTTAAACCCTCTTCAGAGCAAGCGTGGAGGATGATGGAGAATCGTGTGATCAGTGTGCTAAATCTCTCTGCCTGTTTGGACTTTGTAATTATTTTTTTAGCAGTAATTAAAGAAAAAAGTCCTCTGTGAGGAATATTCTCTATTTTAAATATTTTTAGTATGTACTGTGTATGATTCATTACCATTTTGAGGGGATTTATACATATTTTTAGATAAAATTAAATGCTCTTATTTTTCCAACAGCTAAACTACTCTTAGTTGAACAGTGTGCCCTAGCTTTTCTTGCAACCAGAGTATTTTTGTACAGATTTGCTTTCTCTTACAAAAAGAAAAAAAAAATCCTGTTGTATTAACATTTAAAAACAGAATTGTGTTATGTGATCAGTTTTGGGGGTTAACTTTGCTTAATTCCTCAGGCTTTGCGATTTAAGGAGGAGCTGCCTTAAAAAAAAATAAAGGCCTTATTTTGCAATTATGGGAGTAAACAATAGTCTAGAGAAGCATTTGGTAAGCTTTATCATATATATATTTTTTAAAGAAGAGAAAAACACCTTGAGCCTTAAAACGGTGCTGCTGGGAAACATTTGCACTCTTTTAGTGCATTTCCTCCTGCCTTTGCTTGTTCACTGCAGTCTTAAGAAAGAGGTAAAAGGCAAGCAAAGGAGATGAAATCTGTTCTGGGAATGTTTCAGCAGCCAATAAGTGCCCGAGCACACTGCCCCCGGTTGCCTGCCTGGGCCCCATGTGGAAGGCAGATGCCTGCTCGCTCTGTCACCTGTGCCTCTCAGAACACCAGCAGTTAACCTTCAAGACATTCCACTTGCTAAAATTATTTATTTTGTAAGGAGAGGTTTTAATTAAAACAAAAAAAAATTCTTTTTTTTTTTTTTTTCCAATTTTACCTTCTTTAAAATAGGTTGTTGGAGCTTTCCTCAAAGGGTATGGTCATCTGTTGTTAAATTATGTTCTTAACTGTAACCAGTTTTTTTTTATTTATCTCTTTAATCTTTTTTTATTATTAAAAGCAAGTTTCTTTGTATTCCTCACCCTAGATTTGTATAAATGCCTTTTTGTCCATCCCTTTTTTCTTTGTTGTTTTTGTTGAAAACAAACTGGAAACTTGTTTCTTTTTTTGTATAAATGAGAGATTGCAAATGTAGTGTATCACTGAGTCATTTGCAGTGTTTTCTGCCACAGACCTTTGGGCTGCCTTATATTGTGTGTGTGTGTGGGTGTGTGTGTGTTTTGACACAAAAACAATGCAAGCATGTGTCATCCATATTTCTCTGCATCTTCTCTTGGAGTGAGGGAGGCTACCTGGAGGGGATCAGCCCACTGACAGACCTTAATCTTAATTACTGCTGTGGCTAGAGAGTTTGAGGATTGCTTTTTAAAAAAGACAGCAAACTTTTTTTTTTATTTAAAAAAAGATATATTAACAGTTTTAGAAGTCAGTAGAATAAAATCTTAAAGCACTCATAATATGGCATCCTTCAATTTCTGTATAAAAGCAGATCTTTTTAAAAAGATACTTCTGTAACTTAAGAAACCTGGCATTTAAATCATATTTTGTCTTTAGGTAAAAGCTTTGGTTTGTGTTCGTGTTTTGTTTGTTTCACTTGTTTCCCTCCCAGCCCCAAACCTTTTGTTCTCTCCGTGAAACTTACCTTTCCCTTTTTCTTTCTCTTTTTTTTTTTTGTATATTATTGTTTACAATAAATATACATTGCATTAAAAAGAA") 

    def test_get_utrs_negative_strand(self):
        transcript = Transcript("TOR1A", 0)
        both_utr = transcript.get_utr("both")
        print both_utr
        both_utr = transcript.get_seq_from_pos(both_utr)
        self.assertEquals(both_utr, "GAACCGGAAGCGTGGGTCTGGCGGCTGCACCGGTTCGCGGTCGGCGCGAGAACAAGCAGGGTGGCGCGGGTCCGGGCCAGTCATGATTGGCAGCCGGAGTCACTGCCTGGAGTTGGAAAAGAAACAACACTCAGTCCTTCCACACTTCCACCCCCAGCTCCTTTCCCTGGAAGAGGAATCCAGTGAATGTTCCTGTTTGATGTGACAGGAATTCTCCCTGGCATTGTTTCCACCCCCTGGTGCCTGCAGGCCACCCAGGGACCACGGGCGAGGACGTGAAGCCTCCCGAACACGCACAGAAGGAAGGAGCCAGCTCCCAGCCCACTCATCGCAGGGCTCATGATTTTTTACAAATTATGTTTTAATTCCAAGTGTTTCTGTTTCAAGGAAGGATGAATAAGTTTTATTGAAAATGTGGTAACTTTATTTAAAATGATTTTTAACATTATGAGAGACTGCTCAGATTCTAAGTTGTTGGCCTTGTGTGTGTGTTTTTTTTTAAGTTCTCATCATTATTACATAGACTGTGATGTATCTTTACTGGAAATGAGCCCAAGCACACATGCATGGCATTTGTTCCACAGGAGGGCATCCCTGGGGATGTGGCTGGAGCATGAGCCAGCTCTGTCCCAGGATGGTCCCAGCGGATGCTGCCAGGGGCAGTGAAGTGTTTAGGTGAAGGACAAGTAGGTAAGAGGACGCCTTCAGGCACCACAGATAAGCCTGAAACAGCCTCTCCAAGGGTTTTCACCTTAGCAACAATGGGAGCTGTGGGAGTGATTTTGGCCACACTGTCAACATTTGTTAGAACCAGTCTTTTGAAAGAAAAGTATTTCCAACTTGTCACTTGCCAGTCACTCCGTTTTGCAAAAGGTGGCCCTTCACTGTCCATTCCAAATAGCCCACACGTGCTCTCTGCTGGATTCTAAATTATGTGAATTTTGCCATATTAAATCTTCCTCATTTATACTATTATTTGTTACGTTCAATCAGAATCCCCGAAACCTCCTATAAAGCTTAGCTGCCCCTTCTGAGGATGCTGAGAACGGTGTCTTTCTTTATAAATGCAAATGGCTACCGTTTTACAATAAAATTTTGCATGTGCCA") 
        prime_5 = transcript.get_utr("5_prime")
        prime_5 = transcript.get_seq_from_pos(prime_5)
        self.assertEquals(prime_5,"GAACCGGAAGCGTGGGTCTGGCGGCTGCACCGGTTCGCGGTCGGCGCGAGAACAAGCAGGGTGGCGCGGGTCCGGGC") 
        prime_3 = transcript.get_utr("3_prime")
        prime_3 = transcript.get_seq_from_pos(prime_3)
        self.assertEquals(prime_3,"CAGTCATGATTGGCAGCCGGAGTCACTGCCTGGAGTTGGAAAAGAAACAACACTCAGTCCTTCCACACTTCCACCCCCAGCTCCTTTCCCTGGAAGAGGAATCCAGTGAATGTTCCTGTTTGATGTGACAGGAATTCTCCCTGGCATTGTTTCCACCCCCTGGTGCCTGCAGGCCACCCAGGGACCACGGGCGAGGACGTGAAGCCTCCCGAACACGCACAGAAGGAAGGAGCCAGCTCCCAGCCCACTCATCGCAGGGCTCATGATTTTTTACAAATTATGTTTTAATTCCAAGTGTTTCTGTTTCAAGGAAGGATGAATAAGTTTTATTGAAAATGTGGTAACTTTATTTAAAATGATTTTTAACATTATGAGAGACTGCTCAGATTCTAAGTTGTTGGCCTTGTGTGTGTGTTTTTTTTTAAGTTCTCATCATTATTACATAGACTGTGATGTATCTTTACTGGAAATGAGCCCAAGCACACATGCATGGCATTTGTTCCACAGGAGGGCATCCCTGGGGATGTGGCTGGAGCATGAGCCAGCTCTGTCCCAGGATGGTCCCAGCGGATGCTGCCAGGGGCAGTGAAGTGTTTAGGTGAAGGACAAGTAGGTAAGAGGACGCCTTCAGGCACCACAGATAAGCCTGAAACAGCCTCTCCAAGGGTTTTCACCTTAGCAACAATGGGAGCTGTGGGAGTGATTTTGGCCACACTGTCAACATTTGTTAGAACCAGTCTTTTGAAAGAAAAGTATTTCCAACTTGTCACTTGCCAGTCACTCCGTTTTGCAAAAGGTGGCCCTTCACTGTCCATTCCAAATAGCCCACACGTGCTCTCTGCTGGATTCTAAATTATGTGAATTTTGCCATATTAAATCTTCCTCATTTATACTATTATTTGTTACGTTCAATCAGAATCCCCGAAACCTCCTATAAAGCTTAGCTGCCCCTTCTGAGGATGCTGAGAACGGTGTCTTTCTTTATAAATGCAAATGGCTACCGTTTTACAATAAAATTTTGCATGTGCCA") 

if __name__ == '__main__':
    unittest.main()