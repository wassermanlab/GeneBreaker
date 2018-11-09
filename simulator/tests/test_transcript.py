# python -m unittest tests.test_gene 
import unittest
from simulator.transcript import Transcript
from Bio.Seq import Seq

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

    # def test_get_some_info(self):
    #     TOR1A = Transcript("TOR1A", 0)
    #     print "========================tor1a======================="
    #     print "utr 5"
    #     print TOR1A.get_utr("5_prime")
    #     print "coding exons"
    #     print TOR1A.get_coding()
    #     print "introns"
    #     print TOR1A.get_introns()
    #     print "utr 3"
    #     print TOR1A.get_utr("3_prime")        
        
    #     SOX9 = Transcript("SOX9", 0)
    #     print "========================sox9========================"
    #     print "utr 5"
    #     print SOX9.get_utr("5_prime")
    #     print "coding exons"
    #     print SOX9.get_coding()
    #     print "introns"
    #     print SOX9.get_introns()
    #     print "utr 3"
    #     print SOX9.get_utr("3_prime")

    def test_get_sequence_stranded(self):
        TOR1A = Transcript("TOR1A", 0)
        sequence_positive = TOR1A.get_seq()
        sequence_negative = TOR1A.get_seq(True)
        self.assertEquals(Seq(sequence_positive).reverse_complement(), sequence_negative)

    def test_get_codon_positive_easy(self):
        sox9 = Transcript("SOX9", 0)
        codon = sox9.get_codon_from_pos(70117532)
        self.assertEquals(codon[0], "ATG")
        self.assertEquals(codon[1], 0)
        self.assertEquals(codon[2], "+")
    
    def test_get_codon_positive_medium(self):
        sox9 = Transcript("SOX9", 0)
        codon = sox9.get_codon_from_pos(70118859L)
        self.assertEquals(codon[0], "AGA")
        self.assertEquals(codon[1], 2)
        self.assertEquals(codon[2], "+")
    
    def test_get_codon_negative_easy(self):
        sox18 = Transcript("SOX18", 0)
        codon = sox18.get_codon_from_pos(62680869)
        self.assertEquals(codon[0], "ATG")
        self.assertEquals(codon[1], 0)
        self.assertEquals(codon[2], "-")
    
    def test_get_codon_negative_hard(self):
        sox18 = Transcript("SOX18", 0)
        codon = sox18.get_codon_from_pos(62680315)
        self.assertEquals(codon[0], "GGC")
        self.assertEquals(codon[1], 1)
        self.assertEquals(codon[2], "-")


if __name__ == '__main__':
    unittest.main()