# python -m unittest tests.test_gene
import unittest
from MenDelSIM.src.short_tandem_repeat import ShortTandemRepeat
from MenDelSIM.src.transcript import Transcript


class STRBasicTests(unittest.TestCase):
    # test 1
    def test_correct_motif(self):
        STR = {
            "IMPACT": {
                "CHROM": "chr2",
                "END": 191745646,
                "START": 191745598,
                "STR": 5
            },
            "REGION": "UTR",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        STR = ShortTandemRepeat(STR)
        self.assertEqual(STR.get_str_motif(), "GCA")

    def test_retraction_too_large(self):
        STR = {
            "IMPACT": {
                "CHROM": "chr2",
                "END": 191745646,
                "START": 191745598,
                "STR": -50
            },
            "REGION": "GENIC",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        STR = ShortTandemRepeat(STR)
        with self.assertRaises(Exception) as cm:
            STR.get_retraction("chr2")
        err = cm.exception
        self.assertEqual(
            str(err), 'retraction length is larger than the total str')

    def test_retraction_works(self):
        STR = {
            "IMPACT": {
                "CHROM": "chr2",
                "END": 191745646,
                "START": 191745598,
                "STR": -5
            },
            "REGION": "UTR",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        STR = ShortTandemRepeat(STR)
        retraction = STR.get_retraction("chr2")
        self.assertEqual(retraction['pos'], 191745597)
        self.assertEqual(retraction['ref'], "CGCAGCAGCAGCAGCA")
        self.assertEqual(retraction['alt'], "C")

    def test_insertion_too_large(self):
        STR = {
            "IMPACT": {
                "CHROM": "chr2",
                "END": 191745646,
                "START": 191745598,
                "STR": 10000
            },
            "REGION": "UTR",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        STR = ShortTandemRepeat(STR)
        with self.assertRaises(Exception) as cm:
            STR.get_expantion()
        err = cm.exception
        self.assertEqual(str(err), 'expansion length is too large')

    def test_insertion_works(self):
        STR = {
            "IMPACT": {
                "CHROM": "chr2",
                "END": 191745646,
                "START": 191745598,
                "STR": 5
            },
            "REGION": "UTR",
            "TYPE": "STR",
            "ZYGOSITY": "HETEROZYGOUS"}
        STR = ShortTandemRepeat(STR)
        expansion = STR.get_expantion()
        self.assertEqual(expansion['pos'], 191745598)
        self.assertEqual(expansion['ref'], "G")
        self.assertEqual(expansion['alt'], "GCAGCAGCAGCAGCAG")


if __name__ == '__main__':
    unittest.main()