import unittest
from MenDelSIM.src.copy_number_variant import CopyNumberVariant
from MenDelSIM.src.transcript import Transcript
from MenDelSIM.src.api_helper import *


class CNVCreationTests(unittest.TestCase):
    XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
    transcript = Transcript(XKR8_uid, "hg38")

    # test 1
    def test_copy_number_error(self):
        cnv = {
            "TYPE": "CNV",
            "REGION": "GENIC",
            "IMPACT": {
                "START": 27959462,
                "END": 27959865,
                "COPY_CHANGE": 0,
            },
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            CopyNumberVariant(cnv, self.transcript)
        cnv["IMPACT"]["COPY_CHANGE"] = 0
        with self.assertRaises(ValueError):
            CopyNumberVariant(cnv, self.transcript)
        cnv["IMPACT"]["COPY_CHANGE"] = -3
        with self.assertRaises(ValueError):
            CopyNumberVariant(cnv, self.transcript)

    # test 2
    def test_coding_overlap(self):
        cnv = {
            "TYPE": "CNV",
            "REGION": "CODING",
            "IMPACT": {
                "START": 27960050,
                "END": 27960150,
                "COPY_CHANGE": 2,
            },
            "ZYGOSITY": "HETEROZYGOUS"}
        CopyNumberVariant(cnv, self.transcript)

    # test 4
    def test_UTR_overlap(self):
        cnv = {
            "TYPE": "CNV",
            "REGION": "UTR",
            "IMPACT": {
                "START": 27959622,
                "END": 27960150,
                "COPY_CHANGE": 2,
            },
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            CopyNumberVariant(cnv, self.transcript)
        cnv["IMPACT"]["END"] = 27959682
        CopyNumberVariant(cnv, self.transcript)

    # test 5
    def test_intron_overlap(self):
        cnv = {
            "TYPE": "CNV",
            "REGION": "INTRONIC",
            "IMPACT": {
                "START": 27959999,
                "END": 27961050,
                "COPY_CHANGE": 2,
            },
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            CopyNumberVariant(cnv, self.transcript)
        cnv["IMPACT"]["START"] = 27961000
        CopyNumberVariant(cnv, self.transcript)

    # test 6
    def test_deletion(self):
        cnv = {
            "TYPE": "CNV",
            "REGION": "INTRONIC",
            "IMPACT": {
                "START": 27961000,
                "END": 27961050,
                "COPY_CHANGE": -1,
            },
            "ZYGOSITY": "HETEROZYGOUS"}
        cnv = CopyNumberVariant(cnv, self.transcript)
        row = cnv.get_vcf_row()
        self.assertEqual(row['chrom'], "chr1")
        self.assertEqual(row['pos'], "27960999")
        self.assertEqual(
            row['ref'], "GCAACCTCAGCTGGCTTCTTGACCTGGGCCTCCCTGGGAGTTTCAGCAGCCC")
        self.assertEqual(row['alt'], "G")

    # test 7
    def test_insertion(self):
        cnv = {
            "TYPE": "CNV",
            "REGION": "INTRONIC",
            "IMPACT": {
                "START": 27961000,
                "END": 27961050,
                "COPY_CHANGE": 2,
            },
            "ZYGOSITY": "HETEROZYGOUS"}
        cnv = CopyNumberVariant(cnv, self.transcript)
        row = cnv.get_vcf_row()
        self.assertEqual(row['chrom'], "chr1")
        self.assertEqual(row['pos'], "27961000")
        self.assertEqual(
            row['ref'], "CAACCTCAGCTGGCTTCTTGACCTGGGCCTCCCTGGGAGTTTCAGCAGCCC")
        self.assertEqual(
            row['alt'], "CAACCTCAGCTGGCTTCTTGACCTGGGCCTCCCTGGGAGTTTCAGCAGCCCCAACCTCAGCTGGCTTCTTGACCTGGGCCTCCCTGGGAGTTTCAGCAGCCC")
