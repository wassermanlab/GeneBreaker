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

    simple_gene = Gene("simulator/tests/testing_data/genes/basic_gene.json")
    # test 10
    def test_basic_non_coding(self):
        snv = {
        "TYPE": "SNV",
        "REGION": "INTRONIC",
        "IMPACT": {"TYPE_IMPACT": "A", "LOCATION": "ANY"}}
        row = SNV(snv).get_vcf_row(self.simple_gene)
        row = row.split("\t")
        self.assertEqual(row[4].rstrip(), "A")
    # test 11
    def test_basic_non_coding_exact(self):
        snv = {
        "TYPE": "SNV",
        "REGION": "INTRONIC",
        "IMPACT": {"TYPE_IMPACT": "A", "LOCATION": 870}}
        row = SNV(snv).get_non_coding_SNV(self.simple_gene)
        self.assertEqual(row["pos"], 870)
        self.assertEqual(row["ref"], "C")
        self.assertEqual(row["alt"], "A")

class DirectedSNVCodingRegionTests(unittest.TestCase):
    simple_gene = Gene("simulator/tests/testing_data/genes/basic_gene.json")
    # test 12
    def test_directed_SNV_exists_single_replacement(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"TYPE_IMPACT": "G", "LOCATION": 0}}
        snv = SNV(snv)
        self.assertEquals(snv.get_directed_coding_SNV(self.simple_gene, 0)['ref'], "A")
        self.assertEquals(snv.get_directed_coding_SNV(self.simple_gene, 0)['alt'], "G")
    # test 13
    def test_directed_SNV_exists_nonsense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"TYPE_IMPACT": "NONSENSE", "LOCATION": 6}}
        snv = SNV(snv)
        self.assertEquals(snv.get_directed_coding_SNV(self.simple_gene, 6)['ref'], "C")
        self.assertEquals(snv.get_directed_coding_SNV(self.simple_gene, 6)['alt'], "T")
    # test 14
    def test_directed_SNV_exists_missense(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"TYPE_IMPACT": "MISSENSE", "LOCATION": 0}}
        snv = SNV(snv)
        self.assertNotEqual(snv.get_directed_coding_SNV(self.simple_gene, 0), False)
        print snv.get_directed_coding_SNV(self.simple_gene, 0)
    # test 15
    def test_directed_SNV_exists_silent(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"TYPE_IMPACT": "SILENT", "LOCATION": 5}}
        snv = SNV(snv)
        self.assertEquals(snv.get_directed_coding_SNV(self.simple_gene, 5)['ref'], "T")
        self.assertIn(snv.get_directed_coding_SNV(self.simple_gene, 5)['alt'], ["C", "A", "G"])
    # test 16
    def test_directed_SNV_false(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"TYPE_IMPACT": "NONSENSE", "LOCATION": 0}}
        snv = SNV(snv)
        self.assertFalse(snv.get_directed_coding_SNV(self.simple_gene, 0))

class UndirectedSNVCodingRegionTests(unittest.TestCase):
    simple_gene = Gene("simulator/tests/testing_data/genes/basic_gene.json")
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    amino_acid_codons = {"Phe": ["TTT", "TTC"], "Tyr": ["TAT", "TAC"], 
    "His": ["CAT", "CAC"], "Gln": ["CAA", "CAG"], "Asn": ["AAT", "AAC"], 
    "Lys": ["AAA", "AAG"], "Asp": ["GAT", "GAC"], "Glu": ["GAA", "GAG"], 
    "Cys": ["TGT", "TGC"], "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"], 
    "Leu": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], 
    "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], 
    "Pro": ["CCT", "CCC", "CCA", "CCG"], 
    "Gly": ["GGT", "GGC", "GGA", "GGG"], 
    "Val": ["GTT", "GTC", "GTA", "GTG"], 
    "Thr": ["ACT", "ACC", "ACA", "ACG"], 
    "Ala": ["GCT", "GCC", "GCA", "GCG"], 
    "Ile": ["ATT", "ATC", "ATA"], "Trp": ["TGG"]}
    aa_codes = {}
    for key, value in amino_acid_codons.iteritems():
        for code in value:
            aa_codes[code] = key
    # test 17
    def test_undirected_silent_SNV(self):
        snv = {
            "TYPE": "SNV",
            "REGION": "CODING",
            "IMPACT": {"TYPE_IMPACT": "SILENT", "LOCATION": "ANY"}}
        snv = SNV(snv)
        res = snv.get_random_coding_SNV(self.simple_gene)
        self.assertNotIn(res["pos"], range(96,8333))
        codon = self.simple_gene.get_codon_from_pos(res["pos"])
        alt_codon = list(codon[0])
        alt_codon[codon[1]] = res["alt"]
        alt_codon = "".join(alt_codon)
        codon_name = self.aa_codes[codon[0]]
        self.assertIn(alt_codon, self.amino_acid_codons[codon_name])
