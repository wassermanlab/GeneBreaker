# python -m unittest tests.test_gene
import unittest
from MenDelSIM.src.mei import MEI
from MenDelSIM.src.transcript import Transcript
from . import establish_GUD_session
from GUD.ORM import Gene

class MEICreationTests(unittest.TestCase):
    # test 1
    def test_type_value_MEI(self):
        mei = {"TYPE": "SNV",
                 "REGION": "INTRONIC",
                 "IMPACT": {"ELEMENT": "ALU_MELT", "LOCATION": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            MEI(mei)
    # test 2
    def test_location_0(self):
        session = establish_GUD_session()
        SOX18_uid = Gene().select_by_name(session, "SOX18", True)[0].qualifiers["uid"]
        SOX18 = Transcript(SOX18_uid)  # SOX18
        mei = {"TYPE": "MEI",
                 "REGION": "INTRONIC",
                 "IMPACT": {"ELEMENT": "ALU_MELT", "LOCATION": 0},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            MEI(mei).get_vcf_row(SOX18)

class InsersionMethodTesting(unittest.TestCase):
    mei_any = {'TYPE': 'MEI',
                 'REGION': 'INTRONIC',
                 'IMPACT': {"ELEMENT": "ALU_MELT", "LOCATION": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
    mei_any = MEI(mei_any)
    mei_spec = {'TYPE': 'MEI',
                  'REGION': 'CODING',
                  'IMPACT': {"ELEMENT": "ALU_MELT", "LOCATION": 132576250},
                  "ZYGOSITY": "HETEROZYGOUS"}
    mei_spec = MEI(mei_spec)
    session = establish_GUD_session()
    TOR1A_uid = Gene().select_by_name(session, "TOR1A", True)[0].qualifiers["uid"]
    SOX9_uid = Gene().select_by_name(session, "SOX9", True)[0].qualifiers["uid"]    
    positive_transcript = Transcript(SOX9_uid)  # SOX9
    negative_transcript = Transcript(TOR1A_uid)  # TOR1A

    # test 3
    def test_check_insertion_generation(self):
        self.assertEqual(len(self.mei_any.get_insertion_str()), 281)

    # test 8
    def test_get_insertion_any(self):
        insertion = self.mei_any.get_insertion(self.positive_transcript)
        self.assertEqual(len(insertion["alt"]), 282)

    # test 9
    def test_get_insertion_spec(self):
        insertion = self.mei_spec.get_insertion(self.negative_transcript)
        self.assertEqual(insertion["pos"], 132576250)
        self.assertEqual(len(insertion["alt"]), 282)

    # test 10
    def test_get_insertion_row(self):
        insertion = self.mei_spec.get_vcf_row(self.negative_transcript)
        insertion = insertion.split("\t")
        insertion[4] = insertion[4].rstrip()
        self.assertEqual(insertion[0], "chr9")
        self.assertEqual(insertion[1], '132576251')
        self.assertEqual(len(insertion[4]), 282)
        print(insertion)

if __name__ == '__main__':
    unittest.main()