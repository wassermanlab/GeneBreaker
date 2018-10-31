import unittest
from simulator.single_nucleotide_variant import SingleNucleotideVariant as SNV
from simulator.gene import Gene

class SNVCreationTests(unittest.TestCase):
    # test 1
    def test_wrong_type(self):
        snv = {
        "TYPE": "SNV",
        "REGION": "INTRONIC",
        "IMPACT": {"TYPE_IMPACT": "A", "LOCATION": "ANY"}}
        self.assertRaises(SNV(snv))
    # test 2
    def test_wrong_type_impact(self):
        snv = {
        "TYPE": "SNV",
        "REGION": "INTRONIC",
        "IMPACT": {"TYPE_IMPACT": "WRONG", "LOCATION": "ANY"}}
        self.assertRaises(SNV(snv))

class MutationMethods(unittest.TestCase):
    # test 3
    def test_getting_alternate_codons(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"TYPE_IMPACT": "A", "LOCATION": "ANY"}}
        snv = SNV(snv)
        print snv.get_alternate_codons("ATG", 1)
    # test 4
    def test_nonsense_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"TYPE_IMPACT": "A", "LOCATION": "ANY"}}
        snv = SNV(snv)
        self.assertEquals(snv.nonsense_mutation("TAT", 2), "A")
    # test 5
    def test_nonsense_not_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"TYPE_IMPACT": "A", "LOCATION": "ANY"}}
        snv = SNV(snv)
        self.assertFalse(snv.nonsense_mutation("TAT", 1))
    # test 6
    def test_misssense_false(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"TYPE_IMPACT": "A", "LOCATION": "ANY"}}
        snv = SNV(snv)
        self.assertFalse(snv.missense_mutation("TAT", 2))
    # test 7
    def test_misssense_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"TYPE_IMPACT": "A", "LOCATION": "ANY"}}
        snv = SNV(snv)
        print snv.missense_mutation("CAC", 2)
    # test 8
    def test_silent_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"TYPE_IMPACT": "A", "LOCATION": "ANY"}}
        snv = SNV(snv)
        self.assertFalse(snv.silent_mutation("TGG", 2))
    # test 9
    def test_silent_false(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"TYPE_IMPACT": "A", "LOCATION": "ANY"}}
        snv = SNV(snv)
        self.assertEquals(snv.silent_mutation("CAA", 2), "G")

class NonCodingSNV(unittest.TestCase):
    simple_gene = Gene("tests/testing_data/genes/basic_gene.json")
    # test 6
    def test_basic_non_coding(self):
        snv = {
        "TYPE": "SNV",
        "REGION": "INTRONIC",
        "IMPACT": {"TYPE_IMPACT": "A", "LOCATION": "ANY"}}
        row = SNV(snv).get_vcf_row(self.simple_gene)
        row = row.split("\t")
        self.assertEqual(row[4].rstrip(), "A")
    # test 7 
    def test_basic_non_coding_exact(self):
        snv = {
        "TYPE": "SNV",
        "REGION": "INTRONIC",
        "IMPACT": {"TYPE_IMPACT": "A", "LOCATION": 870}}
        row = SNV(snv).get_non_coding_SNV(self.simple_gene)
        self.assertEqual(row["pos"], 870)
        self.assertEqual(row["ref"], "C")
        self.assertEqual(row["alt"], "A")