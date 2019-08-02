import unittest
from MenDelSIM.src.transcript import Transcript
from Bio.Seq import Seq
from MenDelSIM.src.api_helper import *

class TranscriptCreationTests(unittest.TestCase): 
    # @unittest.skip
    def test_normal_gene(self):
        transcript = Transcript(5, 'hg38')
        print(transcript)
    # @unittest.skip
    def test_assert_raises(self):
        with self.assertRaises(Exception) as cm:
            Transcript(9000000, 'hg38')
        err = cm.exception
        self.assertEqual(str(err), 'cannot make transcript')

class TranscriptMethodTests(unittest.TestCase):
    # @unittest.skip
    def test_get_seq(self):
        XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
        transcript = Transcript(XKR8_uid, "hg38")
        seq = transcript.get_seq()
        self.assertEqual(seq, "ATTCCCTTCCGTATGGAAATGAGGTTAACTCATTTGCAGGAACCTAACGCCCGCCCTCCCCCGACTCGGATCGCTGGCGTGAACCCAGCGGCCCCGTGCTCTCCGGGTATCCCTGGATCGAGGCGCACGCGTTGTTGCAAAGTCTCGCACTTAGGGATTAGAATGTGCCTTTGCGTTTCGTTTTGTTTTGCAACAATACAGGTGTGAGGGATCGCAGCCCTGGTCCAGGCTCTAGTCTCCCCGCCCTCTCGCTTTCTTGGGCCTCCGGGCGCCTCAAAGGTAAAACTAGGCCATTTGAAGCCTGCGAAGCGAAGGGCAGTCGGTCTTCAGGCCTAGCGGGCTTGGCTGGGCTGGGTCGGTGGACCCCTGGGCTGGCGGCGGGAGCAGAGCATCGAGCGGCTTTTAGCTGGGGTTCGGGCGGAGGACCTGGGACAACCGCGGGCCTCGCTCCTCCCCGGACGGGGCGGGACGGGGAAGTCCCGCCCACCAGGCCCAGGCCCCGGGCCGCCCCGAGGGCTGCGCCCACCTCCTTCCTGCCTCGGCAACCCCGGGCCCTGAGGGCAGGCCCCAACCGCGGAGGAGCAGGAGAGGGCGGAGGCCGGCGGGCCATGCCCTGGTCGTCCCGCGGCGCCCTCCTTCGGGACCTGGTCCTGGGCGTGCTGGGCACCGCCGCCTTCCTGCTCGACCTGGGCACCGACCTGTGGGCCGCCGTCCAGTATGCGCTCGGCGGCCGCTACCTGTGGGCGGCGCTGGTGCTGGCGCTGCTGGGCCTGGCCTCCGTGGCGCTGCAGCTCTTCAGCTGGCTCTGGCTGCGCGCTGACCCTGCCGGCCTGCACGGGTCGCAGCCCCCGCGCCGCTGCCTGGCGCTGCTGCATCTCCTGCAGCTGGGTTACCTGTACAGGTGAGTGCTTCGCCCCGGGAGGGGAGGAGTGTCGGAGCCCAGCACCTCCGTCAGCTGGGTCACCTGTACAGGTGACACGCCATATGCCGGGAAGGGGACTGGGAGAGGGAACCAAGTCCTGTGGGCGCCTGGCCTACAGGTGAGCGTGCAGCCCTTTACGGAAAGGGAATCCCGGTAGACTGCCTTCATTCTAGCCCTAGCTTTCCTGGCACAGGCGAAGAAACTGAGGCCCAGGGAAGAAAGGGAGGCGCATGGGCTTTGAAGTCGGACCTCTCAGGGATCCACTTCTGGCTTTGCCCCCTCTCGCCTCACTTCTGCATTTCTCTGAGCCTCAGTTTCCTCCTGCGTAGAAAGGAGGTAGTCACCCTGAGACAGGTGTGACGGGCAGGCCTGTTGTGAAGAGTGAGCGAAAGAATGAGAAAGTGCCTCGCACAGTGCCTGGCCTGTAATGGCAGGGGCTGCGGACTCAGTTCAGCATCCCATGTTTCAGCATGGGCAGGGCCTTGGGAAGTCACTGACCTCTCAGTGGGGCTGTTTGGTTGCTGTTCCCATTGGGATGGGGATAGCGCCTGCCCGTTTGGTGAAGATGCAGTGAAATGACATCTGTCCACTGCCCTGGCTCAGAGTGAGCACTCAGCAACCTCAGCTGGCTTCTTGACCTGGGCCTCCCTGGGAGTTTCAGCAGCCCTGACCCCACCACAGCGCTTGCCCCTCACAAAGGCCGCCTCTCACCCTGTCTGCCTGGTTGTATTCTTCCCCCTGAAACTTGGCTTCCCACCCCAGACCCAGCAAGAGGGCAGAGATGGCCTCTGACTCATCCATGCCCATGGAGCCAGCCAGGGCTGGGCATAGGCTGGGTCCTGGAAAATGTGTGTGGATTCTGCTCCTTTTTGGACCTGGTGGAGCCTCGTAGAGATGGGAAAACAAACCCAAGGCTCAGAGGCATGCGACAAAGTCAGGGCAGAGGGTGGTGTGGGTTGGGATGGGGTGGCCAGGGCCAGGATCTGCTCTTGAATCCGTTTTTACAAACTCTTGAGATCTGGCCTTGACAAAGCCACACGGAGTCTCCAGGAGTGGAAGAGAGGAGCTTTGAATGTGGCAGTTTTGCCTGAGGGATTATAGAGACATCAGAAAACCCCGGGTGGGCATGAAGGGGTGCAGTTGTTAAGTCCTGGGCAGGCCTCCATGAGCCGTGAATCTGTAACTACCTTATGCATCTTTGTGGCCAATATTTCTGCTGAGCTTGAGTCAGACAGACTGTACTAGCTTTGGCTTTGGAAATGGGTCTCTGTTCCTGCCTGTGAATTGGCCGAGGCAGGGGTCCCTGCTAGGTCACTGGGGACCATAGGCAAAGGCTCCTGCTGAGGGAGAACCCCAGCCTGTCCTGCCTGCCTGCTCCTCAGTGTGTTCTTGTGTCCTTCCAGGGGAAGTGACAGTACATGTCTGGTGAATGAGCTTGTGATGAGTGAACGGGGTTCTCTCCCAATAGGCAACTGAGAGCAAAGTTCTGGGTGGAGGGTCCTTGGTGCATCCTCCACTTACTCTTCATGCCAGTGGGTACCCCAGACCTGTCTGTCCAAACTCTGGACTCTCTTAACTCCGCCCCAATCTCCATCCTTAACCACTCAGGGCCAGCCTCTTAAATGGGATCCCTCCTTCCAATCTAACCCTCCTCCAGCCCCTCCCCTGTGTGGCCGGGGCAAGGCAGAGTCTGATTCCATCTGGCCTGTGCCCAGAAGTGCCCAAGAGGCAAATGTCCTCTTGCCCACAGGACCCCAGGAGATCTTGGGCACACTAAAGTTTGAGAACTACTGCCTCAGGCAACTTTATAGCATGCCTTTCCAATTCAGGTTTGGATCCAAATCCTCATTCCAGTACTTATTAGCTGCTTATAATTTTTTGTTTGTTTTTTAAAAAGAGGTAGGGTCTGGTCAGGCTTGGTGGCTCATGCCTGTAATCCTAGCAGTCTGGGAGGCTGACGTGGGTGGATCACTTGAGGTCAGGAGTTTGAGACCAGCCTGCTGAACATGGCAAAAGCTCATCTCTACAAAAATTAGCCAGGTGTGGTGGCGCATGCCTGTAATCCCAGCATCCCAGCTACTCTGGAGGCCGAGGCAGGAGAATCGCTTGAACCTGGGAGGTGGAGGTTGCAGTGAGCCGAGTTTGTGCCACTGCACTCTAGCCTGGGCGACAAGCGACAGATGTAGACTCTGTCTCCAAAAAAAAAGAAAAAAAAAAAGAGAGAGAGAGGCAGGGTCTCTGTTGCCTGGGCTGGAGTGCAGTGGCACAATCAATCATGGCTCACAGTGGCCTTGAACTCCTGGGCTCAAGCTGTCCTCCTGCGTCAGCCTCCTGAGTAGCTTGGACTACAGGCTTGTACCACCACGTCCAACTAATTAAAAAAAATTCTTTGTAGAGATGGGGCTTGCTGTGTTGCCCAGGGTGGTCTCAAACTCCTGGCCTCCAGTGATCCTCTTGCCTTGGCCTGTCAAAATGTTGGGATTACAGGCGCGAGCCACTGCACCCAGACCTTATTAGCTTCTTGATCCTGTTTCCCCATCTGTAAAATGGGCTTTTGTGAGAAGTAAACAGGACAAGGTTTGTGATGAGCTTGGGACAGTGCTTAGCCAGCGGGAGTAGAAGGTGTATTGATCCCCATTTTTGCTTTTTTGTTTTTTTTGAGACGGAGTCTTGTTCTGTCGCCAGGCTGGAGTGCAGTGGCACGATCTCGGCTCACTGTAACCTCCGCCTTGCAGATTCAAGAGATTCTCCTGCCTCAGCCTCCTGAGTAGCTGGGACTACAGGCGCGCACCACCACACCCAGCTAATTTTTTGTATTAGTAGAGACGAGGCTTCACCATGTTAGCCAGGATGGTCTCAATCTCCTGAACTCGTGATCAGCCCGCCTCAGACTCCCAAAGTGCTGGGATTACAGGCATGAGCCACCACACCCGGCCTCGATCCTCATTTTACAGATGAAACAAGTTCAGTGAGGTTAGGACATTGCCGAAGGCCCCATAGCATGGATGTGAGGAGACAGGTTGGAGCCTGTGTCCACTCATTAGATGGGTGGGAGGCTGAGGAATTCACAGGACACTAACCTGGCCCTCTGGGCATTTGTGTGTGGTGCCTAGGTGCGTGCAGGAGCTGCGGCAGGGGCTGCTGGTGTGGCAGCAGGAGGAGCCCTCTGAGTTTGACTTGGCCTACGCCGACTTCCTCGCCCTGGACATCAGCATGCTGCGGCTCTTCGAGACCTTCTTGGAGACGGCACCACAGCTCACGCTGGTGCTGGCCATCATGCTGCAGAGTGGCCGGGCTGAGTACTACCAGTGTGAGTGAAGGCCTGTGGCTGGCCCCCCTGTCGTGGCTTGGTGGGGGGTCTCTCAAATGTTGGAACTGTTTTTAGTCTTTTATAAAGGCTGCTTAGAAAACAGGAGAACAGGCTTTAGTCAGGCAGATCTGGCTTGAAACCTAAAGTCACTCTACAGCTGTTTGAGTTTGGACAAATGCCTTGACCTCTCTGAGTATGTTTGTTCTCATCTGAAATATGGGCTTAAATCCTCCTGCCTCATAAGGTTGATGAAAGGATTAAATGAGGTGATGCAAAGAAAGCCCATTTCCTGGTACATAAGTTCCTGGTACAGAGTCTCACTCTGTCATCGCCCGTAGTGGAGTACAGCAGCCTGATCATGGCTCACTGTAGCCTCCACCTCCCAGGCTCAAGTGATCCTCCTGTCTCAGCCTCCTGAGTAGCTGGGACTACAAGTGTGTATCACCATGCTGGCTAATTTTTCTTTCTTTTTTTTTTTTTTTTTTTTGTAGAGACGAGGTCCCACTTTGTTACCCAGGCTGGTCTTGAACTCCTGAGCTCAAGTAACCTCCTACCTCAGCCTCCTAAGGTGCTGGGATTACAGGTGTCGCCCACTGTGCCTGGCCCACAAGTCTTAATTGTAATTTTTATAATTTTGAAGATAATGAAGTGTGTAAGGTGCCTGATACAGCATAGGTAACTTTTTAGTTAAGAAACAATTTCATACGTTGGGAGGCCAAGGTGCACAGGTCGCTTGAGGCCAGAGTTGGAGACCATCCTGGGTAACATAGTGAGACCTCATCTCTTCAAAATTTTAAAAGAGAAACAAGAAACAATATATTGAATGCCTTCATCCAGTCAGGTTTTCATTGTGCATCTCTTCTGTGCCTGTTACTGTGCTGGGGACACAGCAGTGAACAAGATGAACCCAGCCCCTGCCCTGCCAGGATGATAGACTAAAACAAGTAGCTACGGTTGAGTGACTGGAAAAGGGAGAGGCAGGGAGCAGAGGATCACAGGGCCCCCTAAGCATGGGTGAAGTTTACAGTGAGGAGCTTTGGGAGGGTTTTTGCATGAAAGGAAAGTGACTTGCCCATTTCCACAGACCTGGACAGTGGCCAAGCTAAAGAGGGCCCCCACTCACATCCGACTCAGGGTCCAAGCCTTCCCCTTTGCCTTCCTCCACCGCTGCCATAAATGCCACAGCCTCTCAAGAAACCAGTCCTCACTCTACCGTCACTCGCTGTGTGACCTGGAGAGCCTTCCTGTGGAAGATGGAGGTTGGACTCGATCTCCAAGGGCCCTTTCGTCTTGCTAGTCTGAGTCTATATATTGATTGAAAAAACAATAATAGCGGCTGTCACAATCGAGTGCCAGCTATTAGCCAGGCCCTGTGGGAAGCACTTACAGTCATCATTGCTCATGTTCACAGCAGCCCTGTAGGTTTGTGCTATATTGATCTTCATTTTTAAAGAGGTGCAGAAAGGTGAGTGACTTGCCCTGGGTTACTGGGCACTCACTGGGCACACGTCTTTGTCTGTTGAGGGTTGGGGGTGTCTAGAACCAGGGCCAAGTGCAGACAGTCTGCACTGTGAGTGTGGCAGGGGGTAAGGGGGCGAGACAGATTTTCCCTACTTTTTATTTAGCAAACATCTCTTTCGCTGCTGTTATGTGCCAGGTACTGGGCTGCTGGGGGATCCCAAGTGAGCAGAGTCTGTTTTCTACCCTCGAGGAGCCCAGAGATAAGGAAATAGATAATTACTGTGTGATGAGACTCCAGACAGAGGTGGGTGGCATGTTACAGGGATACCTGACACAGCTGAGTGGGTGGCATGGGGTCAGGGTGGGTGACTCCCGGGGCTCTTCAAGCTGGGGCCTGCGGGGTAATGGAGGAGTAGGCCGGGGAAGTGGGGCCTGTGCTTCTGGCGACCTAACATCCTGGGCAGAGCTTGTCAAGGTCACGGAAGGTGGCATGCTGGGCAAACTGGAAGAACATCAGCCCTGCTGAGACCTGGGCAAGGCTAATGAGCACAGATCTGTTTCCAGGCTAGGAGTGGGTTCTGGGCGCCTAATAGGTTCCCAGTAAACATCTGTCGAATGACTGCGGTGGACGGGGACAGGGGAGGAAGCAGTCGGGAGACGTGCCGGTGCCACTTCCATCCCTCCCTGGAGGCCCTTGTGCTTCTGCTCCAGGATTGGCAACCAGGAGGTGGGATAGGACCCCTGTTGTAAACTTCTGGATCTTGGAAAATCAAGGTTAGGGGTCAGTCCGATGGCCTCTGGTCTTGGAGGGGAGGGGAGGGGAGTGGGGAGACCTTCTGTCCACACAGCAGGGACATCTCTCCTGGCCCTGAGGGCAGCTTTAGCCAGAGCCTGTTGCCACCCCATGGGATTTTCTACTCCCTTTTTGAGACTAAGAAAGACGGAACTGAGACTGTTGATTTTTAGACACAAGGAAGGCTTGTGTTCTTCCCAGAAAGAAAGGGTTATGGGGCCCTGGAAGATGCCCAGGATGCAGCCAAGGCCCAGGGGTGCTGAGTGCTGGGTGGGCCACAGGTCACATCCTCTCATCTTTGACAAATACCACCCGGATAGCACCTACGCAGCAAAGCGCTGTGGTGAGGGGCAGGCTCTATATTAGGAGTCATAGCTTTGTTTTTGGTTGTGGGAACTCACTGGGTCTTCTCTAGCTGCAGATTTGCCTTCAACAAACGAGGGATTCTAATACCTACCCCAGGGTTGCTGGGAGGGCCCCCACGGTACCTGTGACCGCTGGGGAGTGCCAAGCAGGCTGGCCCCAGGACTCACTGTTCCCTCCTACTCTTTGCAGGGGTTGGCATCTGCACATCCTTCCTGGGCATCTCGTGGGCACTGCTCGACTACCACCGGGCCTTGCGCACCTGCCTCCCCTCCAAGCCGCTCCTGGGCCTGGGCTCCTCCGTGATCTACTTCCTGTGGAACCTGCTGCTGCTGTGGCCCCGAGTCCTGGCTGTGGCCCTGTTCTCAGCCCTCTTCCCCAGCTATGTGGCCCTGCACTTCCTGGGCCTGTGGCTGGTACTGCTGCTCTGGGTCTGGCTTCAGGGCACAGACTTCATGCCGGACCCCAGCTCCGAGTGGCTGTACCGGGTGACGGTGGCCACCATCCTCTATTTCTCCTGGTTCAACGTGGCTGAGGGCCGCACCCGAGGCCGGGCCATCATCCACTTCGCCTTCCTCCTGAGTGACAGCATTCTCCTGGTGGCCACCTGGGTGACTCATAGCTCCTGGCTGCCCAGCGGGATTCCACTGCAGCTGTGGCTGCCTGTGGGATGCGGCTGCTTCTTTCTGGGCCTGGCTCTGCGGCTTGTGTACTACCACTGGCTGCACCCTAGCTGCTGCTGGAAGCCCGACCCTGACCAGGTAGACGGGGCCCGGAGTCTGCTTTCTCCAGAGGGGTATCAGCTGCCTCAGAACAGGCGCATGACCCATTTAGCACAGAAGTTTTTCCCCAAGGCTAAGGATGAGGCTGCTTCGCCAGTGAAGGGATAGGTGAACGGCGTCCTTTGAAGCAGGATCAGACCCAGCCAGCAGAGATGGAGAGTGACTCTGTTGGCAGAAGGCAGGCGAGGATAAGCTAACGATGCTGCTGTGGCCTCTATGCACTCAGCAAGAGCGGGACGCCTGTGCTGGGCCGGGCACCAGGGATGGTGCTGAGTCGGGCAGAGGCCTCCTTTCAAGGAGTTCACAGTGAACAAGATGAGAAGGGCTGGGCCCTGGAGGGTCAAGAGCCCCAATTATGTACAAGACACTTTGGGAGGAAAGAAGACTACCTTTTCCCCCTGCCATTGGTATAGCTGGTGCCCCAAAACTTCCACCTCCCTCCCTGGCTACCTCTAAAATGACTGGTATAGGTGCTGCCCCACCCCTTAGCTCCCCTATCCTGGGCTAGGAGGCCACAGGGGCTGTCCTCTAGAATTCTTCCTTCCCTCCCCCACACCATTCATTCAATTCATGAAACAAATCTTTGCCAAGAGCAGTTTATGTGCCAGGAACATCATTCTGTCCTTGCAACCTGGAACAAGACCAGCTACCAGCCTAGCTTCATCCGCTACTTGCACCAACCAGTCCCGGGTTAGATCCCAAATGCTAGAAGCCAGGGATGCCCAACTCTGGGTGGCCCCAGTCAGAACCTCTGGGATCTCAGTGAAGCTGGCCTGGCCTCTGCTCCTGCTCTCAAGGGGCTGCTTTTCAACCAAGAGCCTTGTGAGCCTGGTCTGAGCCTTGCACAGCCACTGAGTATTTTTTTTGCCTTAGCCAGTGTACCTCCTACCTCAGTCTATGTGAGAGGAAGAGAATGTGTGTGCCTGTGGGTCTCTACAAGTGACAGATGTGTTGTTTTCAACAGTATTATTAGGTTATGAATAAAGCCTCATGAAATCCTCCAA".upper())

    # @unittest.skip
    def test_get_exons(self):
        XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
        transcript = Transcript(XKR8_uid, "hg38")
        exons = transcript.get_exons()
        exon_seq = transcript.get_seq_from_pos(exons)
        self.assertEqual(exon_seq, "ATTCCCTTCCGTATGGAAATGAGGTTAACTCATTTGCAGGAACCTAACGCCCGCCCTCCCCCGACTCGGATCGCTGGCGTGAACCCAGCGGCCCCGTGCTCTCCGGGTATCCCTGGATCGAGGCGCACGCGTTGTTGCAAAGTCTCGCACTTAGGGATTAGAATGTGCCTTTGCGTTTCGTTTTGTTTTGCAACAATACAGGTGTGAGGGATCGCAGCCCTGGTCCAGGCTCTAGTCTCCCCGCCCTCTCGCTTTCTTGGGCCTCCGGGCGCCTCAAAGGTAAAACTAGGCCATTTGAAGCCTGCGAAGCGAAGGGCAGTCGGTCTTCAGGCCTAGCGGGCTTGGCTGGGCTGGGTCGGTGGACCCCTGGGCTGGCGGCGGGAGCAGAGCATCGAGCGGCTTTTAGCTGGGGTTCGGGCGGAGGACCTGGGACAACCGCGGGCCTCGCTCCTCCCCGGACGGGGCGGGACGGGGAAGTCCCGCCCACCAGGCCCAGGCCCCGGGCCGCCCCGAGGGCTGCGCCCACCTCCTTCCTGCCTCGGCAACCCCGGGCCCTGAGGGCAGGCCCCAACCGCGGAGGAGCAGGAGAGGGCGGAGGCCGGCGGGCCATGCCCTGGTCGTCCCGCGGCGCCCTCCTTCGGGACCTGGTCCTGGGCGTGCTGGGCACCGCCGCCTTCCTGCTCGACCTGGGCACCGACCTGTGGGCCGCCGTCCAGTATGCGCTCGGCGGCCGCTACCTGTGGGCGGCGCTGGTGCTGGCGCTGCTGGGCCTGGCCTCCGTGGCGCTGCAGCTCTTCAGCTGGCTCTGGCTGCGCGCTGACCCTGCCGGCCTGCACGGGTCGCAGCCCCCGCGCCGCTGCCTGGCGCTGCTGCATCTCCTGCAGCTGGGTTACCTGTACAGGTGCGTGCAGGAGCTGCGGCAGGGGCTGCTGGTGTGGCAGCAGGAGGAGCCCTCTGAGTTTGACTTGGCCTACGCCGACTTCCTCGCCCTGGACATCAGCATGCTGCGGCTCTTCGAGACCTTCTTGGAGACGGCACCACAGCTCACGCTGGTGCTGGCCATCATGCTGCAGAGTGGCCGGGCTGAGTACTACCAGTGGGTTGGCATCTGCACATCCTTCCTGGGCATCTCGTGGGCACTGCTCGACTACCACCGGGCCTTGCGCACCTGCCTCCCCTCCAAGCCGCTCCTGGGCCTGGGCTCCTCCGTGATCTACTTCCTGTGGAACCTGCTGCTGCTGTGGCCCCGAGTCCTGGCTGTGGCCCTGTTCTCAGCCCTCTTCCCCAGCTATGTGGCCCTGCACTTCCTGGGCCTGTGGCTGGTACTGCTGCTCTGGGTCTGGCTTCAGGGCACAGACTTCATGCCGGACCCCAGCTCCGAGTGGCTGTACCGGGTGACGGTGGCCACCATCCTCTATTTCTCCTGGTTCAACGTGGCTGAGGGCCGCACCCGAGGCCGGGCCATCATCCACTTCGCCTTCCTCCTGAGTGACAGCATTCTCCTGGTGGCCACCTGGGTGACTCATAGCTCCTGGCTGCCCAGCGGGATTCCACTGCAGCTGTGGCTGCCTGTGGGATGCGGCTGCTTCTTTCTGGGCCTGGCTCTGCGGCTTGTGTACTACCACTGGCTGCACCCTAGCTGCTGCTGGAAGCCCGACCCTGACCAGGTAGACGGGGCCCGGAGTCTGCTTTCTCCAGAGGGGTATCAGCTGCCTCAGAACAGGCGCATGACCCATTTAGCACAGAAGTTTTTCCCCAAGGCTAAGGATGAGGCTGCTTCGCCAGTGAAGGGATAGGTGAACGGCGTCCTTTGAAGCAGGATCAGACCCAGCCAGCAGAGATGGAGAGTGACTCTGTTGGCAGAAGGCAGGCGAGGATAAGCTAACGATGCTGCTGTGGCCTCTATGCACTCAGCAAGAGCGGGACGCCTGTGCTGGGCCGGGCACCAGGGATGGTGCTGAGTCGGGCAGAGGCCTCCTTTCAAGGAGTTCACAGTGAACAAGATGAGAAGGGCTGGGCCCTGGAGGGTCAAGAGCCCCAATTATGTACAAGACACTTTGGGAGGAAAGAAGACTACCTTTTCCCCCTGCCATTGGTATAGCTGGTGCCCCAAAACTTCCACCTCCCTCCCTGGCTACCTCTAAAATGACTGGTATAGGTGCTGCCCCACCCCTTAGCTCCCCTATCCTGGGCTAGGAGGCCACAGGGGCTGTCCTCTAGAATTCTTCCTTCCCTCCCCCACACCATTCATTCAATTCATGAAACAAATCTTTGCCAAGAGCAGTTTATGTGCCAGGAACATCATTCTGTCCTTGCAACCTGGAACAAGACCAGCTACCAGCCTAGCTTCATCCGCTACTTGCACCAACCAGTCCCGGGTTAGATCCCAAATGCTAGAAGCCAGGGATGCCCAACTCTGGGTGGCCCCAGTCAGAACCTCTGGGATCTCAGTGAAGCTGGCCTGGCCTCTGCTCCTGCTCTCAAGGGGCTGCTTTTCAACCAAGAGCCTTGTGAGCCTGGTCTGAGCCTTGCACAGCCACTGAGTATTTTTTTTGCCTTAGCCAGTGTACCTCCTACCTCAGTCTATGTGAGAGGAAGAGAATGTGTGTGCCTGTGGGTCTCTACAAGTGACAGATGTGTTGTTTTCAACAGTATTATTAGGTTATGAATAAAGCCTCATGAAATCCTCCAA")
    
    # @unittest.skip
    def test_get_coding_exons(self):
        XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
        transcript = Transcript(XKR8_uid, "hg38")
        exons = transcript.get_coding()
        exon_seq = transcript.get_seq_from_pos(exons)
        self.assertEqual(exon_seq, "ATGCCCTGGTCGTCCCGCGGCGCCCTCCTTCGGGACCTGGTCCTGGGCGTGCTGGGCACCGCCGCCTTCCTGCTCGACCTGGGCACCGACCTGTGGGCCGCCGTCCAGTATGCGCTCGGCGGCCGCTACCTGTGGGCGGCGCTGGTGCTGGCGCTGCTGGGCCTGGCCTCCGTGGCGCTGCAGCTCTTCAGCTGGCTCTGGCTGCGCGCTGACCCTGCCGGCCTGCACGGGTCGCAGCCCCCGCGCCGCTGCCTGGCGCTGCTGCATCTCCTGCAGCTGGGTTACCTGTACAGGTGCGTGCAGGAGCTGCGGCAGGGGCTGCTGGTGTGGCAGCAGGAGGAGCCCTCTGAGTTTGACTTGGCCTACGCCGACTTCCTCGCCCTGGACATCAGCATGCTGCGGCTCTTCGAGACCTTCTTGGAGACGGCACCACAGCTCACGCTGGTGCTGGCCATCATGCTGCAGAGTGGCCGGGCTGAGTACTACCAGTGGGTTGGCATCTGCACATCCTTCCTGGGCATCTCGTGGGCACTGCTCGACTACCACCGGGCCTTGCGCACCTGCCTCCCCTCCAAGCCGCTCCTGGGCCTGGGCTCCTCCGTGATCTACTTCCTGTGGAACCTGCTGCTGCTGTGGCCCCGAGTCCTGGCTGTGGCCCTGTTCTCAGCCCTCTTCCCCAGCTATGTGGCCCTGCACTTCCTGGGCCTGTGGCTGGTACTGCTGCTCTGGGTCTGGCTTCAGGGCACAGACTTCATGCCGGACCCCAGCTCCGAGTGGCTGTACCGGGTGACGGTGGCCACCATCCTCTATTTCTCCTGGTTCAACGTGGCTGAGGGCCGCACCCGAGGCCGGGCCATCATCCACTTCGCCTTCCTCCTGAGTGACAGCATTCTCCTGGTGGCCACCTGGGTGACTCATAGCTCCTGGCTGCCCAGCGGGATTCCACTGCAGCTGTGGCTGCCTGTGGGATGCGGCTGCTTCTTTCTGGGCCTGGCTCTGCGGCTTGTGTACTACCACTGGCTGCACCCTAGCTGCTGCTGGAAGCCCGACCCTGACCAGGTAGACGGGGCCCGGAGTCTGCTTTCTCCAGAGGGGTATCAGCTGCCTCAGAACAGGCGCATGACCCATTTAGCACAGAAGTTTTTCCCCAAGGCTAAGGATGAGGCTGCTTCGCCAGTGAAGGGATAG")
    
    @unittest.skip # FIXME: fix this test 
    def test_noncoding_gene(self):
        PRDM15_uid = get_all_transcripts("PRDM15", "hg38")[0]["qualifiers"]["uid"]
        transcript = Transcript(PRDM15_uid, "hg38")
        exons = transcript.get_coding()
        exon_seq = transcript.get_seq_from_pos(exons)
        self.assertEqual(exon_seq, "ATGGCTGAAGATGGGAGCGAAGAGATCATGTTCATCTGGTGTGAAGACTGCAGCCAGTACCACGACTCCGAATGTCCCGAGCTGGGCCCAGTGGTCATGGTCAAAGACTCCTTTGTGTTAAGCAGGGCAAGGTCATCCCTTCCTCCCAACTTGGAGATCAGACGACTGGAAGATGGAGCCGAGGGGGTGTTCGCCATCACTCAGCTCGTCAAGCGGACACAGTTCGGTCCCTTTGAGTCCAGGAGGGTCGCCAAATGGGAAAAGGAGTCTGCATTTCCCCTGAAGGTGTTCCAGAAGGACGGGCACCCCGTGTGCTTCGACACCTCCAACGAGGATGACTGCAACTGGATGATGCTGGTGCGGCCAGCGGCGGAGGCCGAGCACCAGAACCTGACGGCCTACCAGCACGGCAGCGACGTGTACTTCACCACCTCCAGAGACATCCCCCCGGGTACCGAGCTGCGCGTGTGGTATGCGGCCTTCTATGCCAAGAAGATGGACAAGCCCATGCTGAAGCAGGCCGGCTCTGGCGTCCACGCTGCAGGCACCCCAGAAAACAGCGCCCCCGTGGAGTCGGAGCCCAGCCAGTGGGCGTGTAAAGTGTGTTCTGCCACCTTCCTGGAGCTGCAGCTCCTCAATGAACATCTGTTGGGCCACTTAGAACAAGCCAAAAGCCTTCCTCCAGGCAGCCAAAGCGAGGCAGCAGCTCCCGAGAAGGAGCAGGACACACCCCGGGGGGAACCCCCTGCAGTGCCCGAGAGCGAGAATGTTGCCACCAAAGAACAGAAGAAAAAGCCTCGAAGGGGGAGAAAACCCAAAGTGTCCAAAGCTGAGCAGCCTCTAGTCATCGTGGAAGACAAGGAACCCACAGAGCAAGTGGCAGAGATCATTACCGAGGTCCCTCCGGATGAGCCTGTGAGTGCAACGCCAGATGAGCGGATCATGGAGCTGGTTCTGGGGAAGCTGGCCACCACCACCACTGACACCAGCTCGGTTCCAAAGTTCACCCATCATCAGAATAACACCATCACGCTCAAGAGGAGCTTAATTCTCTCAAGCAGACACGGCATCCGGCGCAAGCTCATCAAACAGCTCGGGGAGCACAAGCGGGTTTACCAGTGCAATATCTGCAGCAAGATCTTCCAGAACAGCAGCAACCTGAGCAGGCACGTGCGCTCGCATGGTGACAAGCTGTTTAAGTGCGAAGAGTGTGCAAAATTGTTCAGCCGCAAAGAGAGCCTAAAGCAGCACGTTTCCTACAAGCACAGCAGGAACGAGGTGGACGGCGAGTACAGGTACCGCTGCGGCACTTGTGAGAAGACCTTCCGCATCGAGAGCGCGCTGGAGTTCCACAACTGCAGGACAGATGACAAGACGTTCCAATGTGAGATGTGTTTCAGATTCTTCTCCACCAACAGCAACCTCTCCAAGCACAAGAAGAAGCACGGCGACAAGAAGTTTGCCTGTGAGGTCTGCAGCAAGATGTTCTACCGCAAGGACGTCATGCTGGACCACCAGCGCCGGCACCTGGAAGGAGTGCGGCGAGTGAAGCGAGAGGACCTGGAGGCCGGTGGGGAGAACCTGGTCCGTTACAAGAAGGAGCCTTCCGGGTGCCCGGTGTGTGGCAAGGTGTTCTCCTGCCGGAGCAATATGAACAAGCACCTGCTCACCCACGGCGACAAGAAGTACACCTGCGAGATCTGCGGGCGCAAGTTCTTCCGCGTGGATGTGCTCAGGGACCACATCCATGTCCACTTCAAGGACATCGCGTTGATGGATGACCACCAGAGGGAAGAGTTTATCGGCAAGATCGGGATCTCCTCGGAAGAAAACGATGACAATTCTGACGAGAGCGCAGACTCGGAGCCTCACAAGTACAGCTGCAAGCGGTGCCAGCTCACCTTCGGCCGGGGGAAGGAGTACCTGAAGCACATCATGGAGGTGCACAAGGAGAAGGGCTATGGCTGCAGCATCTGCAACCGGCGCTTTGCACTGAAGGCCACCTACCACGCCCACATGGTCATCCACCGTGAAAACCTGCCGGACCCCAACGTGCAGAAGTACATCCACCCCTGCGAGATCTGCGGGCGGATCTTCAACAGCATCGGGAACCTGGAGCGCCACAAGCTCATCCACACAGGTGTGAAGAGCCACGCCTGCGAGCAGTGTGGGAAGTCCTTTGCCAGGAAGGACATGCTGAAGGAGCACATGCGTGTGCACGACAATGTCCGCGAGTACCTGTGTGCCGAGTGTGGGAAAGGCATGAAGACCAAGCACGCGCTGCGCCACCACATGAAGCTGCACAAGGGCATCAAGGAGTACGAGTGCAAGGAGTGCCACCGCAGGTTCGCGCAGAAGGTCAACATGCTCAAGCACTGCAAGCGGCACACGGGGATTAAAGATTTCATGTGTGAATTGTGTGGGAAGACATTCAGCGAGAGGAACACCATGGAGACCCACAAGCTCATCCACACAGTGGGCAAGCAGTGGACGTGCTCCGTGTGCGACAAGAAGTACGTGACCGAGTACATGCTGCAGAAGCACGTTCAGCTCACACACGACAAGGTGGAGGCGCAGAGCTGCCAGCTGTGCGGGACCAAGGTGTCCACCAGGGCCTCCATGAGCCGACACATGCGGCGCAAGCACCCCGAGGTGCTCGCGGTGAGGATCGATGACCTGGACCACCTCCCGGAGACCACCACCATCGACGCCTCCTCCATTGGCATCGTCCAGCCTGAGCTGACTCTGGAGCAGGAGGATTTGGCCGAAGGGAAGCACGGGAAAGCTGCCAAGCGAAGTCACAAGAGAAAGCAGAAGCCAGAAGAGGAGGCGGGTGCTCCGGTGCCCGAGGACGCCACCTTCAGCGAATACTCAGAGAAAGAGACGGAGTTCACAGGCAGTGTAGGCGACGAGACCAATTCCGCAGTACAGAGCATTCAGCAGGTAGTGGTGACCCTGGGTGACCCAAATGTGACCACACCATCGAGCTCAGTCGGCTTAACCAACATCACCGTGACCCCCATCACCACTGCGGCCGCGACTCAGTTTACCAATCTCCAGCCGGTGGCCGTGGGGCACCTTACCACCCCTGAACGCCAGTTACAGCTGGACAACTCAATCCTGACCGTGACCTTTGATACCGTCAGCGGCTCTGCCATGTTGCACAACCGCCAAAATGACGTCCAGATCCACCCCCAGCCGGAAGCCTCGAACCCACAGTCTGTGGCCCATTTCATCAACCTGACGACCCTGGTCAACTCCATCACGCCCCTGGGGAGCCAGCTTAGTGACCAGCACCCGCTCACGTGGCGGGCAGTGCCCCAGACTGACGTCTTGCCACCCTCGCAGCCGCAGGCACCCCCACAGCAGGCGGCCCAGCCCCAGGTGCAGGCGGAGCAGCAGCAGCAGCAGATGTACAGCTACTGA")
    
    # @unittest.skip
    def test_get_introns(self):
        XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
        transcript = Transcript(XKR8_uid, "hg38")
        introns = transcript.get_introns()
        intron_seq = transcript.get_seq_from_pos(introns)
        self.assertEqual(intron_seq, "GTGAGTGCTTCGCCCCGGGAGGGGAGGAGTGTCGGAGCCCAGCACCTCCGTCAGCTGGGTCACCTGTACAGGTGACACGCCATATGCCGGGAAGGGGACTGGGAGAGGGAACCAAGTCCTGTGGGCGCCTGGCCTACAGGTGAGCGTGCAGCCCTTTACGGAAAGGGAATCCCGGTAGACTGCCTTCATTCTAGCCCTAGCTTTCCTGGCACAGGCGAAGAAACTGAGGCCCAGGGAAGAAAGGGAGGCGCATGGGCTTTGAAGTCGGACCTCTCAGGGATCCACTTCTGGCTTTGCCCCCTCTCGCCTCACTTCTGCATTTCTCTGAGCCTCAGTTTCCTCCTGCGTAGAAAGGAGGTAGTCACCCTGAGACAGGTGTGACGGGCAGGCCTGTTGTGAAGAGTGAGCGAAAGAATGAGAAAGTGCCTCGCACAGTGCCTGGCCTGTAATGGCAGGGGCTGCGGACTCAGTTCAGCATCCCATGTTTCAGCATGGGCAGGGCCTTGGGAAGTCACTGACCTCTCAGTGGGGCTGTTTGGTTGCTGTTCCCATTGGGATGGGGATAGCGCCTGCCCGTTTGGTGAAGATGCAGTGAAATGACATCTGTCCACTGCCCTGGCTCAGAGTGAGCACTCAGCAACCTCAGCTGGCTTCTTGACCTGGGCCTCCCTGGGAGTTTCAGCAGCCCTGACCCCACCACAGCGCTTGCCCCTCACAAAGGCCGCCTCTCACCCTGTCTGCCTGGTTGTATTCTTCCCCCTGAAACTTGGCTTCCCACCCCAGACCCAGCAAGAGGGCAGAGATGGCCTCTGACTCATCCATGCCCATGGAGCCAGCCAGGGCTGGGCATAGGCTGGGTCCTGGAAAATGTGTGTGGATTCTGCTCCTTTTTGGACCTGGTGGAGCCTCGTAGAGATGGGAAAACAAACCCAAGGCTCAGAGGCATGCGACAAAGTCAGGGCAGAGGGTGGTGTGGGTTGGGATGGGGTGGCCAGGGCCAGGATCTGCTCTTGAATCCGTTTTTACAAACTCTTGAGATCTGGCCTTGACAAAGCCACACGGAGTCTCCAGGAGTGGAAGAGAGGAGCTTTGAATGTGGCAGTTTTGCCTGAGGGATTATAGAGACATCAGAAAACCCCGGGTGGGCATGAAGGGGTGCAGTTGTTAAGTCCTGGGCAGGCCTCCATGAGCCGTGAATCTGTAACTACCTTATGCATCTTTGTGGCCAATATTTCTGCTGAGCTTGAGTCAGACAGACTGTACTAGCTTTGGCTTTGGAAATGGGTCTCTGTTCCTGCCTGTGAATTGGCCGAGGCAGGGGTCCCTGCTAGGTCACTGGGGACCATAGGCAAAGGCTCCTGCTGAGGGAGAACCCCAGCCTGTCCTGCCTGCCTGCTCCTCAGTGTGTTCTTGTGTCCTTCCAGGGGAAGTGACAGTACATGTCTGGTGAATGAGCTTGTGATGAGTGAACGGGGTTCTCTCCCAATAGGCAACTGAGAGCAAAGTTCTGGGTGGAGGGTCCTTGGTGCATCCTCCACTTACTCTTCATGCCAGTGGGTACCCCAGACCTGTCTGTCCAAACTCTGGACTCTCTTAACTCCGCCCCAATCTCCATCCTTAACCACTCAGGGCCAGCCTCTTAAATGGGATCCCTCCTTCCAATCTAACCCTCCTCCAGCCCCTCCCCTGTGTGGCCGGGGCAAGGCAGAGTCTGATTCCATCTGGCCTGTGCCCAGAAGTGCCCAAGAGGCAAATGTCCTCTTGCCCACAGGACCCCAGGAGATCTTGGGCACACTAAAGTTTGAGAACTACTGCCTCAGGCAACTTTATAGCATGCCTTTCCAATTCAGGTTTGGATCCAAATCCTCATTCCAGTACTTATTAGCTGCTTATAATTTTTTGTTTGTTTTTTAAAAAGAGGTAGGGTCTGGTCAGGCTTGGTGGCTCATGCCTGTAATCCTAGCAGTCTGGGAGGCTGACGTGGGTGGATCACTTGAGGTCAGGAGTTTGAGACCAGCCTGCTGAACATGGCAAAAGCTCATCTCTACAAAAATTAGCCAGGTGTGGTGGCGCATGCCTGTAATCCCAGCATCCCAGCTACTCTGGAGGCCGAGGCAGGAGAATCGCTTGAACCTGGGAGGTGGAGGTTGCAGTGAGCCGAGTTTGTGCCACTGCACTCTAGCCTGGGCGACAAGCGACAGATGTAGACTCTGTCTCCAAAAAAAAAGAAAAAAAAAAAGAGAGAGAGAGGCAGGGTCTCTGTTGCCTGGGCTGGAGTGCAGTGGCACAATCAATCATGGCTCACAGTGGCCTTGAACTCCTGGGCTCAAGCTGTCCTCCTGCGTCAGCCTCCTGAGTAGCTTGGACTACAGGCTTGTACCACCACGTCCAACTAATTAAAAAAAATTCTTTGTAGAGATGGGGCTTGCTGTGTTGCCCAGGGTGGTCTCAAACTCCTGGCCTCCAGTGATCCTCTTGCCTTGGCCTGTCAAAATGTTGGGATTACAGGCGCGAGCCACTGCACCCAGACCTTATTAGCTTCTTGATCCTGTTTCCCCATCTGTAAAATGGGCTTTTGTGAGAAGTAAACAGGACAAGGTTTGTGATGAGCTTGGGACAGTGCTTAGCCAGCGGGAGTAGAAGGTGTATTGATCCCCATTTTTGCTTTTTTGTTTTTTTTGAGACGGAGTCTTGTTCTGTCGCCAGGCTGGAGTGCAGTGGCACGATCTCGGCTCACTGTAACCTCCGCCTTGCAGATTCAAGAGATTCTCCTGCCTCAGCCTCCTGAGTAGCTGGGACTACAGGCGCGCACCACCACACCCAGCTAATTTTTTGTATTAGTAGAGACGAGGCTTCACCATGTTAGCCAGGATGGTCTCAATCTCCTGAACTCGTGATCAGCCCGCCTCAGACTCCCAAAGTGCTGGGATTACAGGCATGAGCCACCACACCCGGCCTCGATCCTCATTTTACAGATGAAACAAGTTCAGTGAGGTTAGGACATTGCCGAAGGCCCCATAGCATGGATGTGAGGAGACAGGTTGGAGCCTGTGTCCACTCATTAGATGGGTGGGAGGCTGAGGAATTCACAGGACACTAACCTGGCCCTCTGGGCATTTGTGTGTGGTGCCTAGGTGAGTGAAGGCCTGTGGCTGGCCCCCCTGTCGTGGCTTGGTGGGGGGTCTCTCAAATGTTGGAACTGTTTTTAGTCTTTTATAAAGGCTGCTTAGAAAACAGGAGAACAGGCTTTAGTCAGGCAGATCTGGCTTGAAACCTAAAGTCACTCTACAGCTGTTTGAGTTTGGACAAATGCCTTGACCTCTCTGAGTATGTTTGTTCTCATCTGAAATATGGGCTTAAATCCTCCTGCCTCATAAGGTTGATGAAAGGATTAAATGAGGTGATGCAAAGAAAGCCCATTTCCTGGTACATAAGTTCCTGGTACAGAGTCTCACTCTGTCATCGCCCGTAGTGGAGTACAGCAGCCTGATCATGGCTCACTGTAGCCTCCACCTCCCAGGCTCAAGTGATCCTCCTGTCTCAGCCTCCTGAGTAGCTGGGACTACAAGTGTGTATCACCATGCTGGCTAATTTTTCTTTCTTTTTTTTTTTTTTTTTTTTGTAGAGACGAGGTCCCACTTTGTTACCCAGGCTGGTCTTGAACTCCTGAGCTCAAGTAACCTCCTACCTCAGCCTCCTAAGGTGCTGGGATTACAGGTGTCGCCCACTGTGCCTGGCCCACAAGTCTTAATTGTAATTTTTATAATTTTGAAGATAATGAAGTGTGTAAGGTGCCTGATACAGCATAGGTAACTTTTTAGTTAAGAAACAATTTCATACGTTGGGAGGCCAAGGTGCACAGGTCGCTTGAGGCCAGAGTTGGAGACCATCCTGGGTAACATAGTGAGACCTCATCTCTTCAAAATTTTAAAAGAGAAACAAGAAACAATATATTGAATGCCTTCATCCAGTCAGGTTTTCATTGTGCATCTCTTCTGTGCCTGTTACTGTGCTGGGGACACAGCAGTGAACAAGATGAACCCAGCCCCTGCCCTGCCAGGATGATAGACTAAAACAAGTAGCTACGGTTGAGTGACTGGAAAAGGGAGAGGCAGGGAGCAGAGGATCACAGGGCCCCCTAAGCATGGGTGAAGTTTACAGTGAGGAGCTTTGGGAGGGTTTTTGCATGAAAGGAAAGTGACTTGCCCATTTCCACAGACCTGGACAGTGGCCAAGCTAAAGAGGGCCCCCACTCACATCCGACTCAGGGTCCAAGCCTTCCCCTTTGCCTTCCTCCACCGCTGCCATAAATGCCACAGCCTCTCAAGAAACCAGTCCTCACTCTACCGTCACTCGCTGTGTGACCTGGAGAGCCTTCCTGTGGAAGATGGAGGTTGGACTCGATCTCCAAGGGCCCTTTCGTCTTGCTAGTCTGAGTCTATATATTGATTGAAAAAACAATAATAGCGGCTGTCACAATCGAGTGCCAGCTATTAGCCAGGCCCTGTGGGAAGCACTTACAGTCATCATTGCTCATGTTCACAGCAGCCCTGTAGGTTTGTGCTATATTGATCTTCATTTTTAAAGAGGTGCAGAAAGGTGAGTGACTTGCCCTGGGTTACTGGGCACTCACTGGGCACACGTCTTTGTCTGTTGAGGGTTGGGGGTGTCTAGAACCAGGGCCAAGTGCAGACAGTCTGCACTGTGAGTGTGGCAGGGGGTAAGGGGGCGAGACAGATTTTCCCTACTTTTTATTTAGCAAACATCTCTTTCGCTGCTGTTATGTGCCAGGTACTGGGCTGCTGGGGGATCCCAAGTGAGCAGAGTCTGTTTTCTACCCTCGAGGAGCCCAGAGATAAGGAAATAGATAATTACTGTGTGATGAGACTCCAGACAGAGGTGGGTGGCATGTTACAGGGATACCTGACACAGCTGAGTGGGTGGCATGGGGTCAGGGTGGGTGACTCCCGGGGCTCTTCAAGCTGGGGCCTGCGGGGTAATGGAGGAGTAGGCCGGGGAAGTGGGGCCTGTGCTTCTGGCGACCTAACATCCTGGGCAGAGCTTGTCAAGGTCACGGAAGGTGGCATGCTGGGCAAACTGGAAGAACATCAGCCCTGCTGAGACCTGGGCAAGGCTAATGAGCACAGATCTGTTTCCAGGCTAGGAGTGGGTTCTGGGCGCCTAATAGGTTCCCAGTAAACATCTGTCGAATGACTGCGGTGGACGGGGACAGGGGAGGAAGCAGTCGGGAGACGTGCCGGTGCCACTTCCATCCCTCCCTGGAGGCCCTTGTGCTTCTGCTCCAGGATTGGCAACCAGGAGGTGGGATAGGACCCCTGTTGTAAACTTCTGGATCTTGGAAAATCAAGGTTAGGGGTCAGTCCGATGGCCTCTGGTCTTGGAGGGGAGGGGAGGGGAGTGGGGAGACCTTCTGTCCACACAGCAGGGACATCTCTCCTGGCCCTGAGGGCAGCTTTAGCCAGAGCCTGTTGCCACCCCATGGGATTTTCTACTCCCTTTTTGAGACTAAGAAAGACGGAACTGAGACTGTTGATTTTTAGACACAAGGAAGGCTTGTGTTCTTCCCAGAAAGAAAGGGTTATGGGGCCCTGGAAGATGCCCAGGATGCAGCCAAGGCCCAGGGGTGCTGAGTGCTGGGTGGGCCACAGGTCACATCCTCTCATCTTTGACAAATACCACCCGGATAGCACCTACGCAGCAAAGCGCTGTGGTGAGGGGCAGGCTCTATATTAGGAGTCATAGCTTTGTTTTTGGTTGTGGGAACTCACTGGGTCTTCTCTAGCTGCAGATTTGCCTTCAACAAACGAGGGATTCTAATACCTACCCCAGGGTTGCTGGGAGGGCCCCCACGGTACCTGTGACCGCTGGGGAGTGCCAAGCAGGCTGGCCCCAGGACTCACTGTTCCCTCCTACTCTTTGCAG".upper())
    
    # @unittest.skip
    def test_get_utrs_positive_strand(self):
        XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
        transcript = Transcript(XKR8_uid, "hg38")
        both_utr = transcript.get_utr("both")
        both_utr = transcript.get_seq_from_pos(both_utr)
        self.assertEqual(both_utr, "ATTCCCTTCCGTATGGAAATGAGGTTAACTCATTTGCAGGAACCTAACGCCCGCCCTCCCCCGACTCGGATCGCTGGCGTGAACCCAGCGGCCCCGTGCTCTCCGGGTATCCCTGGATCGAGGCGCACGCGTTGTTGCAAAGTCTCGCACTTAGGGATTAGAATGTGCCTTTGCGTTTCGTTTTGTTTTGCAACAATACAGGTGTGAGGGATCGCAGCCCTGGTCCAGGCTCTAGTCTCCCCGCCCTCTCGCTTTCTTGGGCCTCCGGGCGCCTCAAAGGTAAAACTAGGCCATTTGAAGCCTGCGAAGCGAAGGGCAGTCGGTCTTCAGGCCTAGCGGGCTTGGCTGGGCTGGGTCGGTGGACCCCTGGGCTGGCGGCGGGAGCAGAGCATCGAGCGGCTTTTAGCTGGGGTTCGGGCGGAGGACCTGGGACAACCGCGGGCCTCGCTCCTCCCCGGACGGGGCGGGACGGGGAAGTCCCGCCCACCAGGCCCAGGCCCCGGGCCGCCCCGAGGGCTGCGCCCACCTCCTTCCTGCCTCGGCAACCCCGGGCCCTGAGGGCAGGCCCCAACCGCGGAGGAGCAGGAGAGGGCGGAGGCCGGCGGGCCGTGAACGGCGTCCTTTGAAGCAGGATCAGACCCAGCCAGCAGAGATGGAGAGTGACTCTGTTGGCAGAAGGCAGGCGAGGATAAGCTAACGATGCTGCTGTGGCCTCTATGCACTCAGCAAGAGCGGGACGCCTGTGCTGGGCCGGGCACCAGGGATGGTGCTGAGTCGGGCAGAGGCCTCCTTTCAAGGAGTTCACAGTGAACAAGATGAGAAGGGCTGGGCCCTGGAGGGTCAAGAGCCCCAATTATGTACAAGACACTTTGGGAGGAAAGAAGACTACCTTTTCCCCCTGCCATTGGTATAGCTGGTGCCCCAAAACTTCCACCTCCCTCCCTGGCTACCTCTAAAATGACTGGTATAGGTGCTGCCCCACCCCTTAGCTCCCCTATCCTGGGCTAGGAGGCCACAGGGGCTGTCCTCTAGAATTCTTCCTTCCCTCCCCCACACCATTCATTCAATTCATGAAACAAATCTTTGCCAAGAGCAGTTTATGTGCCAGGAACATCATTCTGTCCTTGCAACCTGGAACAAGACCAGCTACCAGCCTAGCTTCATCCGCTACTTGCACCAACCAGTCCCGGGTTAGATCCCAAATGCTAGAAGCCAGGGATGCCCAACTCTGGGTGGCCCCAGTCAGAACCTCTGGGATCTCAGTGAAGCTGGCCTGGCCTCTGCTCCTGCTCTCAAGGGGCTGCTTTTCAACCAAGAGCCTTGTGAGCCTGGTCTGAGCCTTGCACAGCCACTGAGTATTTTTTTTGCCTTAGCCAGTGTACCTCCTACCTCAGTCTATGTGAGAGGAAGAGAATGTGTGTGCCTGTGGGTCTCTACAAGTGACAGATGTGTTGTTTTCAACAGTATTATTAGGTTATGAATAAAGCCTCATGAAATCCTCCAA") 
        prime_5 = transcript.get_utr("5_prime")
        prime_5 = transcript.get_seq_from_pos(prime_5)
        self.assertEqual(prime_5,"ATTCCCTTCCGTATGGAAATGAGGTTAACTCATTTGCAGGAACCTAACGCCCGCCCTCCCCCGACTCGGATCGCTGGCGTGAACCCAGCGGCCCCGTGCTCTCCGGGTATCCCTGGATCGAGGCGCACGCGTTGTTGCAAAGTCTCGCACTTAGGGATTAGAATGTGCCTTTGCGTTTCGTTTTGTTTTGCAACAATACAGGTGTGAGGGATCGCAGCCCTGGTCCAGGCTCTAGTCTCCCCGCCCTCTCGCTTTCTTGGGCCTCCGGGCGCCTCAAAGGTAAAACTAGGCCATTTGAAGCCTGCGAAGCGAAGGGCAGTCGGTCTTCAGGCCTAGCGGGCTTGGCTGGGCTGGGTCGGTGGACCCCTGGGCTGGCGGCGGGAGCAGAGCATCGAGCGGCTTTTAGCTGGGGTTCGGGCGGAGGACCTGGGACAACCGCGGGCCTCGCTCCTCCCCGGACGGGGCGGGACGGGGAAGTCCCGCCCACCAGGCCCAGGCCCCGGGCCGCCCCGAGGGCTGCGCCCACCTCCTTCCTGCCTCGGCAACCCCGGGCCCTGAGGGCAGGCCCCAACCGCGGAGGAGCAGGAGAGGGCGGAGGCCGGCGGGCC") 
        prime_3 = transcript.get_utr("3_prime")
        prime_3 = transcript.get_seq_from_pos(prime_3)
        self.assertEqual(prime_3,"GTGAACGGCGTCCTTTGAAGCAGGATCAGACCCAGCCAGCAGAGATGGAGAGTGACTCTGTTGGCAGAAGGCAGGCGAGGATAAGCTAACGATGCTGCTGTGGCCTCTATGCACTCAGCAAGAGCGGGACGCCTGTGCTGGGCCGGGCACCAGGGATGGTGCTGAGTCGGGCAGAGGCCTCCTTTCAAGGAGTTCACAGTGAACAAGATGAGAAGGGCTGGGCCCTGGAGGGTCAAGAGCCCCAATTATGTACAAGACACTTTGGGAGGAAAGAAGACTACCTTTTCCCCCTGCCATTGGTATAGCTGGTGCCCCAAAACTTCCACCTCCCTCCCTGGCTACCTCTAAAATGACTGGTATAGGTGCTGCCCCACCCCTTAGCTCCCCTATCCTGGGCTAGGAGGCCACAGGGGCTGTCCTCTAGAATTCTTCCTTCCCTCCCCCACACCATTCATTCAATTCATGAAACAAATCTTTGCCAAGAGCAGTTTATGTGCCAGGAACATCATTCTGTCCTTGCAACCTGGAACAAGACCAGCTACCAGCCTAGCTTCATCCGCTACTTGCACCAACCAGTCCCGGGTTAGATCCCAAATGCTAGAAGCCAGGGATGCCCAACTCTGGGTGGCCCCAGTCAGAACCTCTGGGATCTCAGTGAAGCTGGCCTGGCCTCTGCTCCTGCTCTCAAGGGGCTGCTTTTCAACCAAGAGCCTTGTGAGCCTGGTCTGAGCCTTGCACAGCCACTGAGTATTTTTTTTGCCTTAGCCAGTGTACCTCCTACCTCAGTCTATGTGAGAGGAAGAGAATGTGTGTGCCTGTGGGTCTCTACAAGTGACAGATGTGTTGTTTTCAACAGTATTATTAGGTTATGAATAAAGCCTCATGAAATCCTCCAA") 

    # @unittest.skip
    def test_get_codon_positive_easy(self):
        XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
        XKR8 = Transcript(XKR8_uid, "hg38")
        codon = XKR8.get_codon_from_pos(27960069)
        self.assertEqual(codon[0], "ATG")
        self.assertEqual(codon[1], 0)
        self.assertEqual(codon[2], 1)

    # @unittest.skip  
    def test_get_codon_positive_medium(self):
        XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
        XKR8 = Transcript(XKR8_uid, "hg38")
        codon = XKR8.get_codon_from_pos(27963497)
        self.assertEqual(codon[0], "TGC")
        self.assertEqual(codon[1], 0)
        self.assertEqual(codon[2], 1)

    # @unittest.skip    
    def test_get_codon_negative_easy(self):
        SOX18_uid = get_all_transcripts("SOX18", "hg38")[0]["qualifiers"]["uid"]
        SOX18 = Transcript(SOX18_uid, "hg38")

        codon = SOX18.get_codon_from_pos(64049514)
        self.assertEqual(codon[0], "ATG")
        self.assertEqual(codon[1], 2)
        self.assertEqual(codon[2], -1)
        

    # @unittest.skip 
    def test_get_codon_negative_hard(self):
        SOX18_uid = get_all_transcripts("SOX18", "hg38")[0]["qualifiers"]["uid"]
        SOX18 = Transcript(SOX18_uid, "hg38")

        codon = SOX18.get_codon_from_pos(64048169)
        self.assertEqual(codon[0], "GGC")
        self.assertEqual(codon[1], 2)
        self.assertEqual(codon[2], -1)

    # @unittest.skip
    def test_get_exons_negative(self):
        TOR1A_uid = get_all_transcripts("TOR1A", "hg38")[0]["qualifiers"]["uid"]
        transcript = Transcript(TOR1A_uid, "hg38")
        exons = transcript.get_exons()
        exon_seq = transcript.get_seq_from_pos(exons)
        self.assertEqual(exon_seq, "GCACCGGTTCGCGGTCGGCGCGAGAACAAGCAGGGTGGCGCGGGTCCGGGCATGAAGCTGGGCCGGGCCGTGCTGGGCCTGCTGCTGCTGGCGCCGTCCGTGGTGCAGGCGGTGGAGCCCATCAGCCTGGGACTGGCCCTGGCCGGCGTCCTCACCGGCTACATCTACCCGCGTCTCTACTGCCTCTTCGCCGAGTGCTGCGGGCAGAAGCGGAGCCTTAGCCGGGAGGCACTGCAGAAGGATCTGGACGACAACCTCTTTGGACAGCATCTTGCAAAGAAAATCATCTTAAATGCCGTGTTTGGTTTCATAAACAACCCAAAGCCCAAGAAACCTCTCACGCTCTCCCTGCACGGGTGGACAGGCACCGGCAAAAATTTCGTCAGCAAGATCATCGCAGAGAATATTTACGAGGGTGGTCTGAACAGTGACTATGTCCACCTGTTTGTGGCCACATTGCACTTTCCACATGCTTCAAACATCACCTTGTACAAGGATCAGTTACAGTTGTGGATTCGAGGCAACGTGAGTGCCTGTGCGAGGTCCATCTTCATATTTGATGAAATGGATAAGATGCATGCAGGCCTCATAGATGCCATCAAGCCTTTCCTCGACTATTATGACCTGGTGGATGGGGTCTCCTACCAGAAAGCCATGTTCATATTTCTCAGCAATGCTGGAGCAGAAAGGATCACAGATGTGGCTTTGGATTTCTGGAGGAGTGGAAAGCAGAGGGAAGACATCAAGCTCAAAGACATTGAACACGCGTTGTCTGTGTCGGTTTTCAATAACAAGAACAGTGGCTTCTGGCACAGCAGCTTAATTGACCGGAACCTCATTGATTATTTTGTTCCCTTCCTCCCCCTGGAATACAAACACCTAAAAATGTGTATCCGAGTGGAAATGCAGTCCCGAGGCTATGAAATTGATGAAGACATTGTAAGCAGAGTGGCTGAGGAGATGACATTTTTCCCCAAAGAGGAGAGAGTTTTCTCAGATAAAGGCTGCAAAACGGTGTTCACCAAGTTAGATTATTACTACGATGATTGACAGTCATGATTGGCAGCCGGAGTCACTGCCTGGAGTTGGAAAAGAAACAACACTCAGTCCTTCCACACTTCCACCCCCAGCTCCTTTCCCTGGAAGAGGAATCCAGTGAATGTTCCTGTTTGATGTGACAGGAATTCTCCCTGGCATTGTTTCCACCCCCTGGTGCCTGCAGGCCACCCAGGGACCACGGGCGAGGACGTGAAGCCTCCCGAACACGCACAGAAGGAAGGAGCCAGCTCCCAGCCCACTCATCGCAGGGCTCATGATTTTTTACAAATTATGTTTTAATTCCAAGTGTTTCTGTTTCAAGGAAGGATGAATAAGTTTTATTGAAAATGTGGTAACTTTATTTAAAATGATTTTTAACATTATGAGAGACTGCTCAGATTCTAAGTTGTTGGCCTTGTGTGTGTGTTTTTTTTTAAGTTCTCATCATTATTACATAGACTGTGATGTATCTTTACTGGAAATGAGCCCAAGCACACATGCATGGCATTTGTTCCACAGGAGGGCATCCCTGGGGATGTGGCTGGAGCATGAGCCAGCTCTGTCCCAGGATGGTCCCAGCGGATGCTGCCAGGGGCAGTGAAGTGTTTAGGTGAAGGACAAGTAGGTAAGAGGACGCCTTCAGGCACCACAGATAAGCCTGAAACAGCCTCTCCAAGGGTTTTCACCTTAGCAACAATGGGAGCTGTGGGAGTGATTTTGGCCACACTGTCAACATTTGTTAGAACCAGTCTTTTGAAAGAAAAGTATTTCCAACTTGTCACTTGCCAGTCACTCCGTTTTGCAAAAGGTGGCCCTTCACTGTCCATTCCAAATAGCCCACACGTGCTCTCTGCTGGATTCTAAATTATGTGAATTTTGCCATATTAAATCTTCCTCATTTATACTATTATTTGTTACGTTCAATCAGAATCCCCGAAACCTCCTATAAAGCTTAGCTGCCCCTTCTGAGGATGCTGAGAACGGTGTCTTTCTTTATAAATGCAAATGGCTACCGTTTTACAATAAAATTTTGCATGTGCCA")
                                    
    # @unittest.skip
    def test_get_coding_negative_exons(self):
        TOR1A_uid = get_all_transcripts("TOR1A", "hg38")[0]["qualifiers"]["uid"]
        transcript = Transcript(TOR1A_uid, "hg38")
        exons = transcript.get_coding()
        exon_seq = transcript.get_seq_from_pos(exons)
        self.assertEqual(exon_seq, "ATGAAGCTGGGCCGGGCCGTGCTGGGCCTGCTGCTGCTGGCGCCGTCCGTGGTGCAGGCGGTGGAGCCCATCAGCCTGGGACTGGCCCTGGCCGGCGTCCTCACCGGCTACATCTACCCGCGTCTCTACTGCCTCTTCGCCGAGTGCTGCGGGCAGAAGCGGAGCCTTAGCCGGGAGGCACTGCAGAAGGATCTGGACGACAACCTCTTTGGACAGCATCTTGCAAAGAAAATCATCTTAAATGCCGTGTTTGGTTTCATAAACAACCCAAAGCCCAAGAAACCTCTCACGCTCTCCCTGCACGGGTGGACAGGCACCGGCAAAAATTTCGTCAGCAAGATCATCGCAGAGAATATTTACGAGGGTGGTCTGAACAGTGACTATGTCCACCTGTTTGTGGCCACATTGCACTTTCCACATGCTTCAAACATCACCTTGTACAAGGATCAGTTACAGTTGTGGATTCGAGGCAACGTGAGTGCCTGTGCGAGGTCCATCTTCATATTTGATGAAATGGATAAGATGCATGCAGGCCTCATAGATGCCATCAAGCCTTTCCTCGACTATTATGACCTGGTGGATGGGGTCTCCTACCAGAAAGCCATGTTCATATTTCTCAGCAATGCTGGAGCAGAAAGGATCACAGATGTGGCTTTGGATTTCTGGAGGAGTGGAAAGCAGAGGGAAGACATCAAGCTCAAAGACATTGAACACGCGTTGTCTGTGTCGGTTTTCAATAACAAGAACAGTGGCTTCTGGCACAGCAGCTTAATTGACCGGAACCTCATTGATTATTTTGTTCCCTTCCTCCCCCTGGAATACAAACACCTAAAAATGTGTATCCGAGTGGAAATGCAGTCCCGAGGCTATGAAATTGATGAAGACATTGTAAGCAGAGTGGCTGAGGAGATGACATTTTTCCCCAAAGAGGAGAGAGTTTTCTCAGATAAAGGCTGCAAAACGGTGTTCACCAAGTTAGATTATTACTACGATGATTGA")

    # @unittest.skip
    def test_get_utrs_negative_strand(self):
        TOR1A_uid = get_all_transcripts("TOR1A", "hg38")[0]["qualifiers"]["uid"]
        transcript = Transcript(TOR1A_uid, "hg38")
        both_utr = transcript.get_utr("both")
        both_utr = transcript.get_seq_from_pos(both_utr)
        self.assertEqual(both_utr, "GCACCGGTTCGCGGTCGGCGCGAGAACAAGCAGGGTGGCGCGGGTCCGGGCCAGTCATGATTGGCAGCCGGAGTCACTGCCTGGAGTTGGAAAAGAAACAACACTCAGTCCTTCCACACTTCCACCCCCAGCTCCTTTCCCTGGAAGAGGAATCCAGTGAATGTTCCTGTTTGATGTGACAGGAATTCTCCCTGGCATTGTTTCCACCCCCTGGTGCCTGCAGGCCACCCAGGGACCACGGGCGAGGACGTGAAGCCTCCCGAACACGCACAGAAGGAAGGAGCCAGCTCCCAGCCCACTCATCGCAGGGCTCATGATTTTTTACAAATTATGTTTTAATTCCAAGTGTTTCTGTTTCAAGGAAGGATGAATAAGTTTTATTGAAAATGTGGTAACTTTATTTAAAATGATTTTTAACATTATGAGAGACTGCTCAGATTCTAAGTTGTTGGCCTTGTGTGTGTGTTTTTTTTTAAGTTCTCATCATTATTACATAGACTGTGATGTATCTTTACTGGAAATGAGCCCAAGCACACATGCATGGCATTTGTTCCACAGGAGGGCATCCCTGGGGATGTGGCTGGAGCATGAGCCAGCTCTGTCCCAGGATGGTCCCAGCGGATGCTGCCAGGGGCAGTGAAGTGTTTAGGTGAAGGACAAGTAGGTAAGAGGACGCCTTCAGGCACCACAGATAAGCCTGAAACAGCCTCTCCAAGGGTTTTCACCTTAGCAACAATGGGAGCTGTGGGAGTGATTTTGGCCACACTGTCAACATTTGTTAGAACCAGTCTTTTGAAAGAAAAGTATTTCCAACTTGTCACTTGCCAGTCACTCCGTTTTGCAAAAGGTGGCCCTTCACTGTCCATTCCAAATAGCCCACACGTGCTCTCTGCTGGATTCTAAATTATGTGAATTTTGCCATATTAAATCTTCCTCATTTATACTATTATTTGTTACGTTCAATCAGAATCCCCGAAACCTCCTATAAAGCTTAGCTGCCCCTTCTGAGGATGCTGAGAACGGTGTCTTTCTTTATAAATGCAAATGGCTACCGTTTTACAATAAAATTTTGCATGTGCCA") 
        prime_5 = transcript.get_utr("5_prime")
        prime_5 = transcript.get_seq_from_pos(prime_5)
        self.assertEqual(prime_5,"GCACCGGTTCGCGGTCGGCGCGAGAACAAGCAGGGTGGCGCGGGTCCGGGC") 
        prime_3 = transcript.get_utr("3_prime")
        prime_3 = transcript.get_seq_from_pos(prime_3)
        self.assertEqual(prime_3,"CAGTCATGATTGGCAGCCGGAGTCACTGCCTGGAGTTGGAAAAGAAACAACACTCAGTCCTTCCACACTTCCACCCCCAGCTCCTTTCCCTGGAAGAGGAATCCAGTGAATGTTCCTGTTTGATGTGACAGGAATTCTCCCTGGCATTGTTTCCACCCCCTGGTGCCTGCAGGCCACCCAGGGACCACGGGCGAGGACGTGAAGCCTCCCGAACACGCACAGAAGGAAGGAGCCAGCTCCCAGCCCACTCATCGCAGGGCTCATGATTTTTTACAAATTATGTTTTAATTCCAAGTGTTTCTGTTTCAAGGAAGGATGAATAAGTTTTATTGAAAATGTGGTAACTTTATTTAAAATGATTTTTAACATTATGAGAGACTGCTCAGATTCTAAGTTGTTGGCCTTGTGTGTGTGTTTTTTTTTAAGTTCTCATCATTATTACATAGACTGTGATGTATCTTTACTGGAAATGAGCCCAAGCACACATGCATGGCATTTGTTCCACAGGAGGGCATCCCTGGGGATGTGGCTGGAGCATGAGCCAGCTCTGTCCCAGGATGGTCCCAGCGGATGCTGCCAGGGGCAGTGAAGTGTTTAGGTGAAGGACAAGTAGGTAAGAGGACGCCTTCAGGCACCACAGATAAGCCTGAAACAGCCTCTCCAAGGGTTTTCACCTTAGCAACAATGGGAGCTGTGGGAGTGATTTTGGCCACACTGTCAACATTTGTTAGAACCAGTCTTTTGAAAGAAAAGTATTTCCAACTTGTCACTTGCCAGTCACTCCGTTTTGCAAAAGGTGGCCCTTCACTGTCCATTCCAAATAGCCCACACGTGCTCTCTGCTGGATTCTAAATTATGTGAATTTTGCCATATTAAATCTTCCTCATTTATACTATTATTTGTTACGTTCAATCAGAATCCCCGAAACCTCCTATAAAGCTTAGCTGCCCCTTCTGAGGATGCTGAGAACGGTGTCTTTCTTTATAAATGCAAATGGCTACCGTTTTACAATAAAATTTTGCATGTGCCA") 
    
    # @unittest.skip 
    def test_get_sequence_stranded(self): 
        TOR1A_uid = get_all_transcripts("TOR1A", "hg38")[0]["qualifiers"]["uid"]
        TOR1A = Transcript(TOR1A_uid, "hg38")
        sequence_positive = TOR1A.get_seq()
        sequence_negative = TOR1A.get_seq(True)
        self.assertEqual(Seq(sequence_positive).reverse_complement(), sequence_negative)
    

if __name__ == '__main__':
    unittest.main()