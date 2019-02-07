# python -m unittest tests.test_gene
import unittest
from simulator.mei import MEI
from simulator.transcript import Transcript


class MEICreationTests(unittest.TestCase):
    # test 4
    def test_type_value_MEI(self):
        indel = {"TYPE": "SNV",
                 "REGION": "INTRONIC",
                 "IMPACT": {"ELEMENT": "ALU_MELT", "LOCATION": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            MEI(mei)
    # test 5
    def test_location_0(self):
        sox18 = Transcript(241)  # SOX18
        indel = {"TYPE": "MEI",
                 "REGION": "INTRONIC",
                 "IMPACT": {"ELEMENT": "ALU_MELT", "LOCATION": 0},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            MEI(mei).get_vcf_row(sox18)

class InsersionMethodTesting(unittest.TestCase):
    indel_any = {'TYPE': 'MEI',
                 'REGION': 'INTRONIC',
                 'IMPACT': {"ELEMENT": "ALU_MELT", "LOCATION": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
    indel_any = Indel(indel_any)
    indel_spec = {'TYPE': 'INDEL',
                  'REGION': 'CODING',
                  'IMPACT': {"INDEL_AMOUNT": 8, "LOCATION": 132576250},
                  "ZYGOSITY": "HETEROZYGOUS"}
    indel_spec = Indel(indel_spec)
    positive_transcript = Transcript(64805)  # SOX9
    negative_transcript = Transcript(68960)  # TOR1A

    # test 7
    def test_check_insertion_generation(self):
        self.assertEqual(len(self.indel_any.get_insertion_str(50)), 50)

    # test 8
    def test_get_insertion_any(self):
        insertion = self.indel_any.get_insertion(self.positive_transcript)
        self.assertEqual(len(insertion["alt"]), 6)

    # test 9
    def test_get_insertion_spec(self):
        insertion = self.indel_spec.get_insertion(self.negative_transcript)
        self.assertEqual(insertion["pos"], 132576250)
        self.assertEqual(len(insertion["alt"]), 9)

    # test 10
    def test_get_insertion_row(self):
        insertion = self.indel_spec.get_vcf_row(self.negative_transcript)
        insertion = insertion.split("\t")
        insertion[4] = insertion[4].rstrip()
        self.assertEquals(insertion[0], "chr9")
        self.assertEquals(insertion[1], '132576251')
        self.assertEquals(len(insertion[4]), 9)
        print insertion

if __name__ == '__main__':
    unittest.main()
