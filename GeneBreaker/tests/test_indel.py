# python -m unittest tests.test_gene
import unittest
import time
from GeneBreaker.src.indel import Indel
from GeneBreaker.src.transcript import Transcript
from GeneBreaker.src.api_helper import *

class IndelCreationTests(unittest.TestCase):
    time.sleep(2)
    XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
    transcript = Transcript(XKR8_uid, "hg38")
    SOX18_uid = get_all_transcripts("SOX18", "hg38")[0]["qualifiers"]["uid"]
    SOX18 = Transcript(SOX18_uid, "hg38")
    # test 1
    def test_wrong_keys(self):
        indel = {
            "TYPE": "INDEL",
            "REGION": "INTRONIC",
            "check": 5,
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(KeyError):
            Indel(indel, self.transcript)
    # test 2
    def test_region_value(self):
        indel = {"TYPE": "INDEL",
                 "REGION": "TEST",
                 "IMPACT": {"INDEL_AMOUNT": 5, "START": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Indel(indel, self.transcript)

    # test 3
    def test_type_value(self):
        indel = {"TYPE": "TEST",
                 "REGION": "INTRONIC",
                 "IMPACT": {"INDEL_AMOUNT": 5, "START": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Indel(indel, self.transcript)

    # test 4
    def test_type_value_indel(self):
        indel = {"TYPE": "SNV",
                 "REGION": "INTRONIC",
                 "IMPACT": {"INDEL_AMOUNT": 8332, "START": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Indel(indel, self.transcript)

    # test 5
    def test_location_0(self):
        indel = {"TYPE": "INDEL",
                 "REGION": "INTRONIC",
                 "IMPACT": {"INDEL_AMOUNT": 5, "START": 0},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Indel(indel, self.SOX18).get_vcf_row()

    # test 6
    def test_location_str(self):
        indel = {"TYPE": "INDEL",
                 "REGION": "INTRONIC",
                 "IMPACT": {"INDEL_AMOUNT": 5, "START": "four"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Indel(indel, self.SOX18).get_vcf_row()

class InsersionMethodTesting(unittest.TestCase):
    time.sleep(2)
    indel_any = {'TYPE': 'INDEL',
                 'REGION': 'INTRONIC',
                 'IMPACT': {"INDEL_AMOUNT": 5, "START": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
    indel_spec = {'TYPE': 'INDEL',
                  'REGION': 'CODING',
                  'IMPACT': {"INDEL_AMOUNT": 8, "START": 129814200},
                  "ZYGOSITY": "HETEROZYGOUS"}
    SOX9_uid = get_all_transcripts("SOX9", "hg38")[0]["qualifiers"]["uid"]
    positive_transcript = Transcript(SOX9_uid, "hg38")
    TOR1A_uid = get_all_transcripts("TOR1A", "hg38")[0]["qualifiers"]["uid"]
    negative_transcript = Transcript(TOR1A_uid, "hg38")
    indel_any = Indel(indel_any, positive_transcript)
    indel_spec = Indel(indel_spec, negative_transcript)
    

    # test 7
    def test_check_insertion_generation(self):
        self.assertEqual(len(self.indel_any.get_insertion_str(50)), 50)

    # test 8
    def test_get_insertion_any(self):
        insertion = self.indel_any.get_insertion()
        self.assertEqual(len(insertion["alt"]), 6)

    # test 9
    def test_get_insertion_spec(self):
        insertion = self.indel_spec.get_insertion()
        self.assertEqual(insertion["pos"], 129814199)
        self.assertEqual(len(insertion["alt"]), 9)

    # test 10
    def test_get_insertion_row(self): 
        insertion = self.indel_spec.get_vcf_row()
        self.assertEqual(insertion['chrom'], "chr9")
        self.assertEqual(insertion['pos'], '129814200')
        self.assertEqual(len(insertion['alt']), 9)
        print(insertion)


class DeletionMethodTestCase(unittest.TestCase):
    time.sleep(2)
    # positive_transcript = Transcript(64805)  # SOX9
    SOX9_uid = get_all_transcripts("SOX9", "hg38")[0]["qualifiers"]["uid"]
    positive_transcript = Transcript(SOX9_uid, "hg38")  
    # test 11
    def test_basic_deletion(self):
        indel = {'TYPE': 'INDEL',
                 'REGION': 'CODING',
                 'IMPACT': {"INDEL_AMOUNT": -5, "START": 72121500},
                 "ZYGOSITY": "HETEROZYGOUS"}
        indel = Indel(indel, self.positive_transcript)
        deletion = indel.get_deletion()
        self.assertEqual(deletion["ref"], "TCGGGC")
        self.assertEqual(deletion["alt"], "T")
        self.assertEqual(deletion["pos"], 72121499)

    # test 12
    def test_any_location_deletion(self):
        indel = {'TYPE': 'INDEL',
                 'REGION': 'CODING',
                 'IMPACT': {"INDEL_AMOUNT": -5, "START": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        indel = Indel(indel, self.positive_transcript)
        deletion = indel.get_deletion()
        self.assertEqual(len(deletion["ref"]), 6)
        self.assertEqual(len(deletion["alt"]), 1)
    
    # test 13
    def test_raises_length_over_200(self):
        indel = {'TYPE': 'INDEL',
                 'REGION': 'CODING',
                 'IMPACT': {"INDEL_AMOUNT": 570, "START": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Indel(indel, self.positive_transcript).get_deletion()

    # test 14
    def test_raises_length_greater_than_region(self): ##########################################
        SOX18_uid = get_all_transcripts("SOX18", "hg38")[0]["qualifiers"]["uid"]
        SOX18 = Transcript(SOX18_uid, "hg38") 
        indel = {'TYPE': 'INDEL',
                 'REGION': 'INTRONIC',
                 'IMPACT': {"INDEL_AMOUNT": -199, "START": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Indel(indel, SOX18)

    # test 15
    def test_raises_length_greater_than_location(self):
        SOX18_uid = get_all_transcripts("SOX18", "hg38")[0]["qualifiers"]["uid"]
        SOX18 = Transcript(SOX18_uid, "hg38")
        indel = {'TYPE': 'INDEL',
                 'REGION': 'INTRONIC',
                 'IMPACT': {"INDEL_AMOUNT": -150, "START": 64049150},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            indel = Indel(indel, SOX18)
            


if __name__ == '__main__':
    unittest.main()
