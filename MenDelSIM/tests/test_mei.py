# python -m unittest tests.test_gene
import unittest
from MenDelSIM.src.mei import MEI
from MenDelSIM.src.transcript import Transcript
from MenDelSIM.src.api_helper import *


class MEICreationTests(unittest.TestCase):
    # test 1
    def test_type_value_MEI(self):
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
        SOX18_uid = get_all_transcripts("SOX18", "hg38")[0]["qualifiers"]["uid"]
        SOX18 = Transcript(SOX18_uid, "hg38")
        mei = {"TYPE": "MEI",
                 "REGION": "INTRONIC",
                 "IMPACT": {"ELEMENT": "ALU_MELT", "START": 0},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            MEI(mei,SOX18).get_vcf_row()

class InsersionMethodTesting(unittest.TestCase):
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

    # test 3
    def test_check_insertion_generation(self):
        mei = MEI(self.mei_any, self.positive_transcript)
        self.assertEqual(len(mei.get_insertion_str()), 281)

    # test 4
    def test_get_insertion_any(self):
        mei = MEI(self.mei_any, self.positive_transcript)
        insertion = mei.get_insertion()
        self.assertEqual(len(insertion["alt"]), 282)

    # test 5
    def test_get_insertion_spec(self):
        mei = MEI(self.mei_spec, self.negative_transcript)
        insertion = mei.get_insertion()
        self.assertEqual(insertion["pos"], 129813999)
        self.assertEqual(len(insertion["alt"]), 282)

    # test 6
    def test_get_insertion_row(self):
        mei = MEI(self.mei_spec, self.negative_transcript)
        insertion = mei.get_vcf_row()
        insertion = insertion.split("\t")
        insertion[4] = insertion[4].rstrip()
        self.assertEqual(insertion[0], "chr9")
        self.assertEqual(insertion[1], '129814000')
        self.assertEqual(len(insertion[4]), 282)
        print(insertion)

if __name__ == '__main__':
    unittest.main()
