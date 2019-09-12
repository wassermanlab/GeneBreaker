import unittest
from MenDelSIM.src.single_nucleotide_variant import SingleNucleotideVariant as SNV
from MenDelSIM.src.transcript import Transcript
from MenDelSIM.src.api_helper import *


class SNVCreationTests(unittest.TestCase):
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
        print(self.snv.missense_mutation("CAC", 2))

    # test 7
    def test_synonymous_exists(self):
        self.assertFalse(self.snv.synonymous_mutation("TGG", 2))

    # test 8
    def test_synonymous_false(self):
        self.assertEqual(self.snv.synonymous_mutation("CAA", 2), "G")

    # test 9
    def test_any_mutation(self):
        self.assertIn(self.snv.any_mutation("CAA", 2), ["G", "C", "T"])


class NonCodingSNV(unittest.TestCase):
    SOX9_uid = get_all_transcripts("SOX9", "hg38")[0]["qualifiers"]["uid"]
    positive_transcript = Transcript(SOX9_uid, "hg38")
    # test 10

    def test_basic_non_coding(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "START": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.positive_transcript)
        row = snv.get_vcf_row()
        self.assertEqual(row["alt"], "A")

    # test 11
    def test_basic_non_coding_exact(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "START": 72122000},
            "ZYGOSITY": "HETEROZYGOUS"}
        row = SNV(snv, self.positive_transcript).get_non_coding_SNV()
        self.assertEqual(row["pos"], 72121999)
        self.assertEqual(row["ref"], "C")
        self.assertEqual(row["alt"], "A")

    # test 12
    def test_non_coding_incorrect_impact(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "MISSENSE", "START": 72121999},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            snv = SNV(snv, self.positive_transcript)
            snv.get_non_coding_SNV()


class DirectedSNVCodingRegionTestsPositive(unittest.TestCase):
    SOX9_uid = get_all_transcripts("SOX9", "hg38")[0]["qualifiers"]["uid"]
    positive_transcript = Transcript(SOX9_uid, "hg38")
    # test 13

    def test_directed_SNV_exists_single_replacement(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "G", "START": 72121500},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.positive_transcript)
        directed = snv.get_directed_coding_SNV(72121499)
        self.assertEqual(directed['ref'], "T")
        self.assertEqual(directed['alt'], "G")

    # test 14
    def test_directed_SNV_exists_nonsense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "NONSENSE", "START": 72121513},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.positive_transcript)
        mutation_nonsense = snv.get_directed_coding_SNV(72121512)
        self.assertEqual(mutation_nonsense['ref'], "C")
        self.assertEqual(mutation_nonsense['alt'], "A")

    # test 15
    def test_directed_SNV_exists_missense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "MISSENSE", "START": 72121517},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.positive_transcript)
        directed_missense = snv.get_directed_coding_SNV(72121516)
        self.assertNotEqual(directed_missense, False)
        print(directed_missense)

    # test 16
    def test_directed_SNV_exists_SYNONYMOUS(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "SYNONYMOUS", "START": 72121517},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.positive_transcript)
        mutation = snv.get_directed_coding_SNV(72121516)
        self.assertEqual(mutation['ref'], "C")
        self.assertEqual(mutation['alt'], "T")

    # test 17
    def test_directed_SNV_false(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "NONSENSE", "START": 72121517},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.positive_transcript)
        self.assertFalse(snv.get_directed_coding_SNV(0))


class DirectedSNVCodingRegionTestsNegative(unittest.TestCase):
    SOX18_uid = get_all_transcripts("SOX18", "hg38")[0]["qualifiers"]["uid"]
    negative_transcript = Transcript(SOX18_uid, "hg38")
    # test 18

    def test_directed_SNV_exists_single_replacement(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "G", "START": 64048520},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.negative_transcript)
        directed = snv.get_directed_coding_SNV(64048519)
        self.assertEqual(directed['ref'], "C")
        self.assertEqual(directed['alt'], "G")

    # test 19
    def test_directed_SNV_exists_nonsense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "NONSENSE", "START": 64048519},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.negative_transcript)
        mutation_nonsense = snv.get_directed_coding_SNV(64048518)
        print(mutation_nonsense)
        self.assertEqual(mutation_nonsense['ref'], "C")
        self.assertEqual(mutation_nonsense['alt'], "A")

    # test 20
    def test_directed_SNV_not_exists_missense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "MISSENSE", "START": 64048502},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.negative_transcript)
        directed_missense = snv.get_directed_coding_SNV(64048501)
        self.assertEqual(directed_missense, False)

    # test 21
    def test_directed_SNV_exists_SYNONYMOUS(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "SYNONYMOUS", "START": 64048502},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv, self.negative_transcript)
        mutation = snv.get_directed_coding_SNV(64048501)
        self.assertEqual(mutation['ref'], "C")
        self.assertIn(mutation['alt'], ["A", "G", "T"])
