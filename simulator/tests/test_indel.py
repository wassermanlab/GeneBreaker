# python -m unittest tests.test_gene
import unittest
from simulator.indel import Indel
from simulator.transcript import Transcript


class IndelCreationTests(unittest.TestCase):
    # test 1
    def test_wrong_keys(self):
        indel = {
            "TYPE": "INDEL",
            "REGION": "INTRONIC",
            "check": 5,
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(KeyError):
            Indel(indel)
    # test 2
    def test_region_value(self):
        indel = {"TYPE": "INDEL",
                 "REGION": "TEST",
                 "IMPACT": {"INDEL_AMOUNT": 5, "LOCATION": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Indel(indel)
    # test 3
    def test_type_value(self):
        indel = {"TYPE": "TEST",
                 "REGION": "INTRONIC",
                 "IMPACT": {"INDEL_AMOUNT": 5, "LOCATION": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Indel(indel)

    # test 4
    def test_type_value_indel(self):
        indel = {"TYPE": "SNV",
                 "REGION": "INTRONIC",
                 "IMPACT": {"INDEL_AMOUNT": 8332, "LOCATION": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Indel(indel)
    # test 5
    def test_location_0(self):
        sox18 = Transcript(241)  # SOX18
        indel = {"TYPE": "INDEL",
                 "REGION": "INTRONIC",
                 "IMPACT": {"INDEL_AMOUNT": 5, "LOCATION": 0},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Indel(indel).get_vcf_row(sox18)
    # test 6
    def test_location_str(self):
        sox18 = Transcript(241)  # SOX18
        indel = {"TYPE": "INDEL",
                 "REGION": "INTRONIC",
                 "IMPACT": {"INDEL_AMOUNT": 5, "LOCATION": "four"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Indel(indel).get_vcf_row(sox18)

class InsersionMethodTesting(unittest.TestCase):
    indel_any = {'TYPE': 'INDEL',
                 'REGION': 'INTRONIC',
                 'IMPACT': {"INDEL_AMOUNT": 5, "LOCATION": "ANY"},
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


class DeletionMethodTestCase(unittest.TestCase):
    positive_transcript = Transcript(64805)  # SOX9
    
    # test 11
    def test_basic_deletion(self):
        indel = {'TYPE': 'INDEL',
                 'REGION': 'CODING',
                 'IMPACT': {"INDEL_AMOUNT": -5, "LOCATION": 70117532},
                 "ZYGOSITY": "HETEROZYGOUS"}
        indel = Indel(indel)
        deletion = indel.get_deletion(self.positive_transcript)
        self.assertEqual(deletion["ref"], "ATGAAT")
        self.assertEqual(deletion["alt"], "A")
        self.assertEqual(deletion["pos"], 70117532)

    # test 12
    def test_any_location_deletion(self):
        indel = {'TYPE': 'INDEL',
                 'REGION': 'CODING',
                 'IMPACT': {"INDEL_AMOUNT": -5, "LOCATION": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        indel = Indel(indel)
        deletion = indel.get_deletion(self.positive_transcript)
        self.assertEqual(len(deletion["ref"]), 6)
        self.assertEqual(len(deletion["alt"]), 1)
    
    # test 13
    def test_raises_length_over_200(self):
        indel = {'TYPE': 'INDEL',
                 'REGION': 'CODING',
                 'IMPACT': {"INDEL_AMOUNT": 570, "LOCATION": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Indel(indel).get_deletion(self.positive_transcript)

    # test 14
    def test_raises_length_greater_than_region(self):
        sox18 = Transcript(241)  # SOX18
        indel = {'TYPE': 'INDEL',
                 'REGION': 'INTRONIC',
                 'IMPACT': {"INDEL_AMOUNT": -199, "LOCATION": "ANY"},
                 "ZYGOSITY": "HETEROZYGOUS"}
        indel = Indel(indel)
        with self.assertRaises(ValueError):
            indel.get_deletion(sox18)

    # test 15
    def test_raises_length_greater_than_location(self):
        sox18 = Transcript(241)  # SOX18
        indel = {'TYPE': 'INDEL',
                 'REGION': 'INTRONIC',
                 'IMPACT': {"INDEL_AMOUNT": -150, "LOCATION": 62680411},
                 "ZYGOSITY": "HETEROZYGOUS"}
        indel = Indel(indel)
        with self.assertRaises(ValueError):
            indel.get_deletion(sox18)


if __name__ == '__main__':
    unittest.main()
