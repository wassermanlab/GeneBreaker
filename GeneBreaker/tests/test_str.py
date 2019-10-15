# python -m unittest tests.test_gene
import unittest
from GeneBreaker.src.short_tandem_repeat import ShortTandemRepeat
from GeneBreaker.src.transcript import Transcript
from GeneBreaker.src.api_helper import *


class STRBasicTests(unittest.TestCase):
    WASH7P_uid = get_all_transcripts("WASH7P", "hg38")[0]["qualifiers"]["uid"]
    transcript = Transcript(WASH7P_uid, "hg38")
    STR_uid = get_strs(16620, 16631, "chr1", "hg38", "exact")[0]["qualifiers"]["uid"]
    # test 1

    def test_correct_motif(self):
        STR = {
            "IMPACT": {
                "STR_ID": self.STR_uid,
                "LENGTH": 2,
            },
            "REGION": "UTR",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        STR = ShortTandemRepeat(STR, self.transcript)
        self.assertEqual(STR.motif, "GCT")

    def test_retraction_too_large(self):
        STR = {
            "IMPACT": {
                "STR_ID": self.STR_uid,
                "LENGTH": -50,
            },
            "REGION": "GENIC",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(Exception) as cm:
            STR = ShortTandemRepeat(STR, self.transcript)
        err = cm.exception
        self.assertEqual(
            str(err), 'retraction length is larger than the total str')

    def test_retraction_works(self):
        STR = {
            "IMPACT": {
                "STR_ID": self.STR_uid,
                "LENGTH": -2,
            },
            "REGION": "UTR",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        STR = ShortTandemRepeat(STR, self.transcript)
        vcf_row = STR.get_vcf_row()
        self.assertEqual(vcf_row['pos'], "16619")
        self.assertEqual(vcf_row['ref'], "CGCTGCT")
        self.assertEqual(vcf_row['alt'], "C")

    def test_insertion_too_large(self):
        STR = {
            "IMPACT": {
                "STR_ID": self.STR_uid,
                "LENGTH": 10000,
            },
            "REGION": "UTR",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(Exception) as cm:
            STR = ShortTandemRepeat(STR, self.transcript)
        err = cm.exception
        self.assertEqual(str(err), 'expansion length is too large')

    def test_insertion_works(self):
        STR = {
            "IMPACT": {
                "STR_ID": self.STR_uid,
                "LENGTH": 5,
            },
            "REGION": "UTR",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        STR = ShortTandemRepeat(STR, self.transcript)
        vcf_row = STR.get_vcf_row()
        self.assertEqual(vcf_row['pos'], "16620")
        self.assertEqual(vcf_row['ref'], "GCT")
        self.assertEqual(vcf_row['alt'], "GCTGCTGCTGCTGCTGCT")


if __name__ == '__main__':
    unittest.main()
