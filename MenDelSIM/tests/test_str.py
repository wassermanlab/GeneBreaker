# python -m unittest tests.test_gene
import unittest
from MenDelSIM.src.short_tandem_repeat import ShortTandemRepeat
from MenDelSIM.src.transcript import Transcript
from MenDelSIM.src.api_helper import *

class STRBasicTests(unittest.TestCase):
    WASH7P_uid = get_all_transcripts("WASH7P", "hg38")[0]["qualifiers"]["uid"]
    transcript = Transcript(WASH7P_uid, "hg38")
    # test 1
    def test_correct_motif(self):
        STR = {
            "IMPACT": {
                "END": 16631,
                "START": 16620,
                "STR": 2, 
                "MOTIF": "GCT"
            },
            "REGION": "UTR",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        STR = ShortTandemRepeat(STR, self.transcript)

    def test_retraction_too_large(self):
        STR = {
            "IMPACT": {
                "END": 16631,
                "START": 16620,
                "STR": -50, 
                "MOTIF": "GCT"
            },
            "REGION": "GENIC",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        STR = ShortTandemRepeat(STR, self.transcript)
        with self.assertRaises(Exception) as cm:
            STR.get_retraction()
        err = cm.exception
        self.assertEqual(
            str(err), 'retraction length is larger than the total str')

    def test_retraction_works(self):
        STR = {
            "IMPACT": {
                "END": 16631,
                "START": 16620,
                "STR": -2, 
                "MOTIF": "GCT"
            },
            "REGION": "UTR",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        STR = ShortTandemRepeat(STR, self.transcript)
        retraction = STR.get_retraction()
        self.assertEqual(retraction['pos'], 16618)
        self.assertEqual(retraction['ref'], "CGCTGCT")
        self.assertEqual(retraction['alt'], "C")

    def test_insertion_too_large(self):
        STR = {
            "IMPACT": {
                "END": 16631,
                "START": 16620,
                "STR": 10000, 
                "MOTIF": "GCT"
            },
            "REGION": "UTR",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        STR = ShortTandemRepeat(STR, self.transcript)
        with self.assertRaises(Exception) as cm:
            STR.get_expantion()
        err = cm.exception
        self.assertEqual(str(err), 'expansion length is too large')

    def test_insertion_works(self):
        STR = {
            "IMPACT": {
                "END": 16631,
                "START": 16620,
                "STR": 5, 
                "MOTIF": "GCT"
            },
            "REGION": "UTR",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        STR = ShortTandemRepeat(STR, self.transcript)
        expansion = STR.get_expantion()
        self.assertEqual(expansion['pos'], 16619)
        self.assertEqual(expansion['ref'], "GCT")
        self.assertEqual(expansion['alt'], "GCTGCTGCTGCTGCTGCT")


if __name__ == '__main__':
    unittest.main()
