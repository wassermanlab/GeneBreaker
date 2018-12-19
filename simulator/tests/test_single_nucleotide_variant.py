import unittest
from simulator.single_nucleotide_variant import SingleNucleotideVariant as SNV
from simulator.transcript import Transcript

class SNVCreationTests(unittest.TestCase):
    # test 1
    def test_wrong_type(self):
        snv = {
            "TYPE": "SNP",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            SNV(snv)
    
    # test 2
    def test_wrong_type_impact(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "WRONG", "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            SNV(snv)

class MutationMethods(unittest.TestCase):
    # test 3
    def test_getting_alternate_codons(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
    
    # test 4
    def test_nonsense_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        self.assertEquals(snv.nonsense_mutation("TAT", 2), "A")
    
    # test 5
    def test_nonsense_not_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        self.assertFalse(snv.nonsense_mutation("TAT", 1))
    
    # test 6
    def test_misssense_false(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        self.assertFalse(snv.missense_mutation("TAT", 2))
    
    # test 7
    def test_misssense_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        print snv.missense_mutation("CAC", 2)
    
    # test 8
    def test_synonymous_exists(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        self.assertFalse(snv.synonymous_mutation("TGG", 2))
    
    # test 9
    def test_synonymous_false(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        self.assertEquals(snv.synonymous_mutation("CAA", 2), "G")
   
    # test 10
    def test_any_mutation(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        self.assertIn(snv.any_mutation("CAA", 2), ["G", "C", "T"])


class NonCodingSNV(unittest.TestCase):
    positive_transcript = Transcript(64805)  # SOX9
    # test 11
    def test_basic_non_coding(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        row = snv.get_vcf_row(self.positive_transcript)
        row = row.split("\t")
        self.assertEqual(row[4].rstrip(), "A")
    
    # test 12
    def test_basic_non_coding_exact(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "A", "LOCATION": 70117963},
            "ZYGOSITY": "HETEROZYGOUS"}
        row = SNV(snv).get_non_coding_SNV(self.positive_transcript)
        self.assertEqual(row["pos"], 70117963)
        self.assertEqual(row["ref"], "G")
        self.assertEqual(row["alt"], "A")
   
    # test 13
    def test_non_coding_incorrect_impact(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "INTRONIC",
            "IMPACT": {"SNV_TYPE": "MISSENSE", "LOCATION": 70117963},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            SNV(snv)
            snv.get_non_coding_SNV(self.positive_transcript)


class DirectedSNVCodingRegionTestsPositive(unittest.TestCase):
    positive_transcript = Transcript(64805)  # SOX9 
    # test 14
    def test_directed_SNV_exists_single_replacement(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "G", "LOCATION": 70117532},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        directed = snv.get_directed_coding_SNV(
            self.positive_transcript, 70117532)
        self.assertEquals(directed['ref'], "A")
        self.assertEquals(directed['alt'], "G")
    
    # test 15
    def test_directed_SNV_exists_nonsense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "NONSENSE", "LOCATION": 70117556},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        mutation_nonsense = snv.get_directed_coding_SNV(
            self.positive_transcript, 70117556)
        self.assertEquals(mutation_nonsense['ref'], "A")
        self.assertEquals(mutation_nonsense['alt'], "T")
    
    # test 16
    def test_directed_SNV_exists_missense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "MISSENSE", "LOCATION": 70117535},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        directed_missense = snv.get_directed_coding_SNV(
            self.positive_transcript, 70117535)
        self.assertNotEqual(directed_missense, False)
        print directed_missense
    
    # test 17
    def test_directed_SNV_exists_SYNONYMOUS(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "SYNONYMOUS", "LOCATION": 70117537},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        mutation = snv.get_directed_coding_SNV(
            self.positive_transcript, 70117537)
        self.assertEquals(mutation['ref'], "T")
        self.assertEquals(mutation['alt'], "C")
    
    # test 18
    def test_directed_SNV_false(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "NONSENSE", "LOCATION": 0},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        self.assertFalse(snv.get_directed_coding_SNV(
            self.positive_transcript, 0))


class DirectedSNVCodingRegionTestsNegative(unittest.TestCase):
    negative_transcript = Transcript(241)  # SOX18
    # test 19
    def test_directed_SNV_exists_single_replacement(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "G", "LOCATION": 62680869},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        directed = snv.get_directed_coding_SNV(
            self.negative_transcript, 62680869)
        self.assertEquals(directed['ref'], "T")
        self.assertEquals(directed['alt'], "G")
    
    # test 20
    def test_directed_SNV_exists_nonsense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "NONSENSE", "LOCATION": 62680866},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        mutation_nonsense = snv.get_directed_coding_SNV(
            self.negative_transcript, 62680866)
        self.assertEquals(mutation_nonsense['ref'], "G")
        self.assertEquals(mutation_nonsense['alt'], "A")
    
    # test 21
    def test_directed_SNV_exists_missense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "MISSENSE", "LOCATION": 62680860},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        directed_missense = snv.get_directed_coding_SNV(
            self.negative_transcript, 62680860)
        self.assertNotEqual(directed_missense, False)
        print directed_missense
    
    # test 22
    def test_directed_SNV_exists_SYNONYMOUS(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "SYNONYMOUS", "LOCATION": 62680846},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        mutation = snv.get_directed_coding_SNV(
            self.negative_transcript, 62680846)
        self.assertEquals(mutation['ref'], "G")
        self.assertEquals(mutation['alt'], "A")
    
    # test 23
    def test_weird(self):
        gene = Transcript(32593)  # TP53
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"SNV_TYPE": "ANY", "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        snv = SNV(snv)
        mutation = snv.get_random_coding_SNV(gene)
        print mutation
