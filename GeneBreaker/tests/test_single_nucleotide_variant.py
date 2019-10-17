import unittest
from GeneBreaker.src.single_nucleotide_variant import SingleNucleotideVariant as SNV
from GeneBreaker.src.transcript import Transcript
from GeneBreaker.src.api_helper import *
import time


class SNVCreationTests(unittest.TestCase):
    time.sleep(1)
    XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
    transcript = Transcript(XKR8_uid, "hg38")
    # test 1

    def test_wrong_type(self):
        snv = {
            "TYPE": "SNP",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "START": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            SNV(snv, self.transcript)

    # test 2
    def test_wrong_type_impact(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "WRONG", "START": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            SNV(snv, self.transcript)


class MutationMethods(unittest.TestCase):
    time.sleep(1)
    XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
    transcript = Transcript(XKR8_uid, "hg38")
    snv = {
        "TYPE": "SNV",
        "REGION": "INTRONIC",
        "IMPACT": {"SNV_TYPE": "A", "START": "ANY"},
        "ZYGOSITY": "HETEROZYGOUS"}
    snv = SNV(snv, transcript)

    # test 3
    def test_nonsense_exists(self):
        self.assertEqual(self.snv.nonsense_mutation("TAT", 2), "A")

    # test 4
    def test_nonsense_not_exists(self):
        self.assertFalse(self.snv.nonsense_mutation("TAT", 1))

    # test 5
    def test_misssense_false(self):
        self.assertFalse(self.snv.missense_mutation("TAT", 2))

    # test 6
    def test_misssense_exists(self):
        self.assertIn(self.snv.missense_mutation("CAC", 2), ["G", "A"])

    # test 7
    def test_synonymous_exists(self):
        self.assertFalse(self.snv.synonymous_mutation("TGG", 2))

    # test 8
    def test_synonymous_false(self):
        self.assertEqual(self.snv.synonymous_mutation("CAA", 2), "G")

    # test 9
    def test_any_mutation(self):
        self.assertIn(self.snv.any_mutation("CAA", 2), ["G", "C", "T"])
    
    # test 10
    def test_stoploss(self):
        self.assertIn(self.snv.stoploss_mutation("TAG", 2), ["C", "T"])
    
    # test 11
    def test_stoploss_not_valid(self):
        self.assertFalse(self.snv.stoploss_mutation("GGA", 2))


class NonCodingSNV(unittest.TestCase):
    time.sleep(1)
    SOX9_uid = get_all_transcripts("SOX9", "hg38")[0]["qualifiers"]["uid"]
    positive_transcript = Transcript(SOX9_uid, "hg38")
    
    # test 12
    def test_basic_non_coding(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "START": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.positive_transcript)
        row = snv.get_vcf_row()
        self.assertEqual(row["alt"], "A")

    # test 13
    def test_basic_non_coding_exact(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "START": 72122000},
            "ZYGOSITY": "HETEROZYGOUS"}
        row = SNV(snv, self.positive_transcript).get_vcf_row()
        self.assertEqual(row["pos"], "72122000")
        self.assertEqual(row["ref"], "C")
        self.assertEqual(row["alt"], "A")

    # test 14
    def test_non_coding_incorrect_impact(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "MISSENSE", "START": 72121999},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            snv = SNV(snv, self.positive_transcript)


class DirectedSNVCodingRegionTestsPositive(unittest.TestCase):
    time.sleep(1)
    SOX9_uid = get_all_transcripts("SOX9", "hg38")[0]["qualifiers"]["uid"]
    positive_transcript = Transcript(SOX9_uid, "hg38")
    
    # test 15
    def test_directed_SNV_exists_single_replacement(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "G", "START": 72121500},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.positive_transcript).get_vcf_row()
        self.assertEqual(snv['ref'], "T")
        self.assertEqual(snv['alt'], "G")

    # test 16
    def test_directed_SNV_exists_nonsense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "NONSENSE", "START": 72121513},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.positive_transcript).get_vcf_row()
        self.assertEqual(snv['ref'], "C")
        self.assertEqual(snv['alt'], "A")

    # test 17
    def test_directed_SNV_exists_missense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "MISSENSE", "START": 72121517},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.positive_transcript)

    # test 18
    def test_directed_SNV_exists_SYNONYMOUS(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "SYNONYMOUS", "START": 72121517},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.positive_transcript).get_vcf_row()
        self.assertEqual(snv['ref'], "C")
        self.assertEqual(snv['alt'], "T")

    # test 19
    def test_directed_SNV_false(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "NONSENSE", "START": 72121517},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(Exception):
            snv = SNV(snv, self.positive_transcript)

    # test 20
    def test_directed_SNV_stoploss_false(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "STOPLOSS", "START": 72124380},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(Exception):
            snv = SNV(snv, self.positive_transcript)
        
    # test 21
    def test_directed_SNV_stoploss_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "STOPLOSS", "START": 72124385},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.positive_transcript).get_vcf_row()
        self.assertEqual(snv['ref'], "T")
        self.assertIn(snv['alt'], ["A","G","C"])
        self.assertEqual(snv['pos'], "72124385")
    
    # test 22
    def test_undirected_SNV_stoploss_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "STOPLOSS", "START": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.positive_transcript).get_vcf_row()
        self.assertIn(snv['pos'], ["72124385", "72124386", "72124387"])



class DirectedSNVCodingRegionTestsNegative(unittest.TestCase):
    time.sleep(1)
    SOX18_uid = get_all_transcripts("SOX18", "hg38")[0]["qualifiers"]["uid"]
    negative_transcript = Transcript(SOX18_uid, "hg38")
    
    # test 23
    def test_directed_SNV_exists_single_replacement(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "G", "START": 64048520},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.negative_transcript).get_vcf_row()
        self.assertEqual(snv['ref'], "C")
        self.assertEqual(snv['alt'], "G")

    # test 24
    def test_directed_SNV_exists_nonsense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "NONSENSE", "START": 64048519},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.negative_transcript).get_vcf_row()
        self.assertEqual(snv['ref'], "C")
        self.assertEqual(snv['alt'], "A")

    # test 25
    def test_directed_SNV_not_exists_missense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "MISSENSE", "START": 64048502},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(Exception):
            snv = SNV(snv, self.negative_transcript)

    # test 26
    def test_directed_SNV_exists_SYNONYMOUS(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "SYNONYMOUS", "START": 64048502},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.negative_transcript).get_vcf_row()
        self.assertEqual(snv['ref'], "C")
        self.assertIn(snv['alt'], ["A", "G", "T"])

    # test 27
    def test_directed_SNV_stoploss_false(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "STOPLOSS", "START": 64048164},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(Exception):
            snv = SNV(snv, self.negative_transcript)
        
    # test 28
    def test_directed_SNV_stoploss_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "STOPLOSS", "START": 64048166},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.negative_transcript).get_vcf_row()
        self.assertEqual(snv['ref'], "C")
        self.assertIn(snv['alt'], ["A","G"])
        self.assertEqual(snv['pos'], "64048166")
    
    # test 29
    def test_undirected_SNV_stoploss_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "STOPLOSS", "START": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.negative_transcript).get_vcf_row()
        self.assertIn(snv['pos'], ["64048166", "64048167", "64048168"])
    
    # test 29
    def test_undirected_SNV(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "ANY", "START": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.negative_transcript).get_vcf_row()

