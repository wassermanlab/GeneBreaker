import unittest
from simulator.single_nucleotide_variant import SingleNucleotideVariant as SNV
from simulator.transcript import Transcript

# @unittest.skip("SNVCreationTests")
class SNVCreationTests(unittest.TestCase):
    # test 1
    def test_wrong_type(self):
        snv = {
        "TYPE": "SNV",
        "REGION": "INTRONIC",
        "IMPACT": "A", "LOCATION": "ANY"}
        self.assertRaises(SNV(snv))
    # test 2
    def test_wrong_type_impact(self):
        snv = {
        "TYPE": "SNV",
        "REGION": "INTRONIC",
        "IMPACT": "WRONG", "LOCATION": "ANY"}
        self.assertRaises(SNV(snv))

# @unittest.skip("MutationMethods")
class MutationMethods(unittest.TestCase):
    # test 3
    def test_getting_alternate_codons(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": "A", "LOCATION": "ANY"}
        snv = SNV(snv)
    # test 4
    def test_nonsense_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": "A", "LOCATION": "ANY"}
        snv = SNV(snv)
        self.assertEquals(snv.nonsense_mutation("TAT", 2), "A")
    # test 5
    def test_nonsense_not_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": "A", "LOCATION": "ANY"}
        snv = SNV(snv)
        self.assertFalse(snv.nonsense_mutation("TAT", 1))
    # test 6
    def test_misssense_false(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": "A", "LOCATION": "ANY"}
        snv = SNV(snv)
        self.assertFalse(snv.missense_mutation("TAT", 2))
    # test 7
    def test_misssense_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": "A", "LOCATION": "ANY"}
        snv = SNV(snv)
        print snv.missense_mutation("CAC", 2)
    # test 8
    def test_silent_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": "A", "LOCATION": "ANY"}
        snv = SNV(snv)
        self.assertFalse(snv.silent_mutation("TGG", 2))
    # test 9
    def test_silent_false(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": "A", "LOCATION": "ANY"}
        snv = SNV(snv)
        self.assertEquals(snv.silent_mutation("CAA", 2), "G")
    # test 10
    def test_any_mutation(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": "A", "LOCATION": "ANY"}
        snv = SNV(snv)
        self.assertIn(snv.any_mutation("CAA", 2), ["G", "C", "T"])

# @unittest.skip("NonCodingSNV")
class NonCodingSNV(unittest.TestCase):
    positive_transcript = Transcript(64805) #SOX9
    # test 11
    def test_basic_non_coding(self):
        snv = {
        "TYPE": "SNV",
        "REGION": "INTRONIC",
        "IMPACT": "A", "LOCATION": "ANY"}
        snv = SNV(snv)
        row = snv.get_vcf_row(self.positive_transcript)
        row = row.split("\t")
        self.assertEqual(row[4].rstrip(), "A")
    # test 12
    def test_basic_non_coding_exact(self):
        snv = {
        "TYPE": "SNV",
        "REGION": "INTRONIC",
        "IMPACT": "A", "LOCATION": 70117963}
        row = SNV(snv).get_non_coding_SNV(self.positive_transcript)
        self.assertEqual(row["pos"], 70117963)
        self.assertEqual(row["ref"], "G")
        self.assertEqual(row["alt"], "A")
    # test 13
    def test_non_coding_incorrect_impact(self):
        snv = {
        "TYPE": "SNV",
        "REGION": "INTRONIC",
        "IMPACT": "MISSENSE", "LOCATION": 70117963}
        snv = SNV(snv)
        print snv.get_non_coding_SNV(self.positive_transcript)
        self.assertRaises(snv.get_non_coding_SNV(self.positive_transcript))

# @unittest.skip("DirectedSNVCodingRegionTestsPositive")
class DirectedSNVCodingRegionTestsPositive(unittest.TestCase):
    positive_transcript = Transcript(64805) #SOX9
    # test 14
    def test_directed_SNV_exists_single_replacement(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": "G", "LOCATION": 70117532}
        snv = SNV(snv)
        directed = snv.get_directed_coding_SNV(self.positive_transcript, 70117532)
        self.assertEquals(directed['ref'], "A")
        self.assertEquals(directed['alt'], "G")
    # test 15
    def test_directed_SNV_exists_nonsense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": "NONSENSE", "LOCATION": 70117556}
        snv = SNV(snv)
        mutation_nonsense = snv.get_directed_coding_SNV(self.positive_transcript, 70117556)
        self.assertEquals(mutation_nonsense['ref'], "A")
        self.assertEquals(mutation_nonsense['alt'], "T")
    # test 16
    def test_directed_SNV_exists_missense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": "MISSENSE", "LOCATION": 70117535}
        snv = SNV(snv)
        directed_missense = snv.get_directed_coding_SNV(self.positive_transcript, 70117535)
        self.assertNotEqual(directed_missense, False)
        print directed_missense
    # test 17
    def test_directed_SNV_exists_silent(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": "SILENT", "LOCATION": 70117537}
        snv = SNV(snv)
        mutation = snv.get_directed_coding_SNV(self.positive_transcript, 70117537)
        self.assertEquals(mutation['ref'], "T")
        self.assertEquals(mutation['alt'], "C")
    # test 18
    def test_directed_SNV_false(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": "NONSENSE", "LOCATION": 0}
        snv = SNV(snv)
        self.assertFalse(snv.get_directed_coding_SNV(self.positive_transcript, 0))

class DirectedSNVCodingRegionTestsNegative(unittest.TestCase):
    negative_transcript = Transcript(241) # SOX18
    # test 19
    def test_directed_SNV_exists_single_replacement(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": "G", "LOCATION": 62680869}
        snv = SNV(snv)
        directed = snv.get_directed_coding_SNV(self.negative_transcript, 62680869)
        self.assertEquals(directed['ref'], "T")
        self.assertEquals(directed['alt'], "G")
    # test 20
    def test_directed_SNV_exists_nonsense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": "NONSENSE", "LOCATION": 62680866}
        snv = SNV(snv)
        mutation_nonsense = snv.get_directed_coding_SNV(self.negative_transcript, 62680866)
        self.assertEquals(mutation_nonsense['ref'], "G")
        self.assertEquals(mutation_nonsense['alt'], "A")
    # test 21
    def test_directed_SNV_exists_missense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": "MISSENSE", "LOCATION": 62680860}
        snv = SNV(snv)
        directed_missense = snv.get_directed_coding_SNV(self.negative_transcript, 62680860)
        self.assertNotEqual(directed_missense, False)
        print directed_missense
    #test 22
    def test_directed_SNV_exists_silent(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": "SILENT", "LOCATION": 62680846}
        snv = SNV(snv)
        mutation = snv.get_directed_coding_SNV(self.negative_transcript, 62680846)
        self.assertEquals(mutation['ref'], "G")
        self.assertEquals(mutation['alt'], "A")
    #test 22
    def test_weird(self):
            gene = Transcript(32593) #TP53
            snv = {
                "TYPE": "SNV",
                "REGION": "CODING",
                "IMPACT": "ANY", "LOCATION": "ANY"}
            snv = SNV(snv)
            mutation = snv.get_random_coding_SNV(gene)
            print mutation
