import unittest
import time 
from GeneBreaker.src.copy_number_variant import CopyNumberVariant
from GeneBreaker.src.transcript import Transcript
from GeneBreaker.src.api_helper import *


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
        with self.assertRaises(ValueError):
            CopyNumberVariant(cnv, self.transcript)

    # test 4
    def test_UTR_overlap(self): ###############################################CHECK####################################
        cnv = {
            "TYPE": "CNV",
            "REGION": "UTR",
            "IMPACT": {
                "START": 27959781,
                "END": 27959832,
                "COPY_CHANGE": 1,
            },
            "ZYGOSITY": "HETEROZYGOUS"}
        cnv["IMPACT"]["END"] = 27959883
        CopyNumberVariant(cnv, self.transcript)

    # test 5
    def test_intron_overlap(self):
        cnv = {
            "TYPE": "CNV",
            "REGION": "INTRONIC",
            "IMPACT": {
                "START": 27959999,
                "END": 27961050,
                "COPY_CHANGE": 1,
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
        cnv = CopyNumberVariant(cnv, self.transcript).get_vcf_row()
        self.assertEqual(cnv['chrom'], "1")
        self.assertEqual(cnv['pos'], "27961000")
        self.assertEqual(cnv['ref'], "C")
        self.assertEqual(cnv['alt'], "<DEL>")

    # test 7
    def test_insertion(self):
        cnv = {
            "TYPE": "CNV",
            "REGION": "INTRONIC",
            "IMPACT": {
                "START": 27961000,
                "END": 27961050,
                "COPY_CHANGE": 1,
            },
            "ZYGOSITY": "HETEROZYGOUS"}
        cnv = CopyNumberVariant(cnv, self.transcript).get_vcf_row()
        self.assertEqual(cnv['chrom'], "1")
        self.assertEqual(cnv['pos'], "27961000")
        self.assertEqual(
            cnv['ref'], "C")
        self.assertEqual(
            cnv['alt'], "<DUP:TANDEM>")
