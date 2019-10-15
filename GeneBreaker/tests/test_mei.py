# python -m unittest tests.test_gene
import unittest
import time
import re
from GeneBreaker.src.mei import MEI
from GeneBreaker.src.transcript import Transcript
from GeneBreaker.src.api_helper import *


class MEICreationTests(unittest.TestCase):
    # test 1
    def test_type_value_MEI(self):
        time.sleep(2)
        XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
        transcript = Transcript(XKR8_uid, "hg38")
        mei = {"TYPE": "SNV",
                 "REGION": "INTRONIC",
                 "IMPACT": {"ELEMENT": "ALU", "START": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            MEI(mei, transcript)
    # test 2
    def test_location_0(self):
        time.sleep(2)
        SOX18_uid = get_all_transcripts("SOX18", "hg38")[0]["qualifiers"]["uid"]
        SOX18 = Transcript(SOX18_uid, "hg38")
        mei = {"TYPE": "MEI",
                 "REGION": "INTRONIC",
                 "IMPACT": {"ELEMENT": "ALU", "START": 0},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            MEI(mei,SOX18).get_vcf_row()

class InsersionMethodTesting(unittest.TestCase):
    time.sleep(2)
    mei_any = {'TYPE': 'MEI',
                 'REGION': 'INTRONIC',
                 'IMPACT': {"ELEMENT": "ALU", "START": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
    mei_spec = {'TYPE': 'MEI',
                  'REGION': 'CODING',
                  'IMPACT': {"ELEMENT": "ALU", "START": 129814000},
                  "ZYGOSITY": "HETEROZYGOUS"}
    SOX9_uid = get_all_transcripts("SOX9", "hg38")[0]["qualifiers"]["uid"]
    positive_transcript = Transcript(SOX9_uid, "hg38")
    TOR1A_uid = get_all_transcripts("TOR1A", "hg38")[0]["qualifiers"]["uid"]
    negative_transcript = Transcript(TOR1A_uid, "hg38")

    # test 4
    def test_get_insertion_any(self):
        mei = MEI(self.mei_any, self.positive_transcript)
        vcf_row = mei.get_vcf_row()
        self.assertEqual(vcf_row["alt"], '<INS:MEI:ALU>')
        self.assertRegex(vcf_row["info"], re.compile('.*SVLEN=282.*'))
    
    # test 5
    def test_get_insertion_spec(self):
        mei = MEI(self.mei_spec, self.negative_transcript)
        vcf_row = mei.get_vcf_row()
        self.assertEqual(vcf_row["pos"], '129814000')
        self.assertEqual(vcf_row["alt"], '<INS:MEI:ALU>')
        self.assertEqual(vcf_row["info"], 'SVTYPE=INS;END=129814282;SVLEN=282;')
        self.assertEqual(vcf_row['chrom'], "chr9")

if __name__ == '__main__':
    unittest.main()
