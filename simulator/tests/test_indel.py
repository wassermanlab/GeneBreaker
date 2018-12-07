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
            "LOCATION": "ANY",
            "ZYGOSITY": "HETEROZYGOUS"}
        self.assertRaises(Indel(indel))
    # test 2
    def test_region_value(self):
        indel = {"TYPE": "INDEL",
                 "REGION": "TEST",
                 "IMPACT": 5,
                 "LOCATION": "ANY",
                 "ZYGOSITY": "HETEROZYGOUS"}
        self.assertRaises(Indel(indel))
    # test 3
    def test_type_value(self):
        indel = {"TYPE": "TEST",
                 "REGION": "INTRONIC",
                 "IMPACT": 5,
                 "LOCATION": "ANY",
                 "ZYGOSITY": "HETEROZYGOUS"}
        self.assertRaises(Indel(indel))
    # test 4
    def test_type_value_indel(self):
        indel = {"TYPE": "SNV",
                 "REGION": "INTRONIC",
                 "IMPACT": 8332,
                 "LOCATION": "ANY"}
        self.assertRaises(Indel(indel))
    # test 5
    def test_location_0(self):
        indel = {"TYPE": "INDEL",
                 "REGION": "INTRONIC",
                 "IMPACT": 5,
                 "LOCATION": 0,
                 "ZYGOSITY": "HETEROZYGOUS"}
        self.assertRaises(Indel(indel))
    # test 6
    def test_location_str(self):
        indel = {"TYPE": "INDEL",
                 "REGION": "INTRONIC",
                 "IMPACT": 5,
                 "LOCATION": "four",
                 "ZYGOSITY": "HETEROZYGOUS"}
        self.assertRaises(Indel(indel))


class InheritanceVariantMethodTests(unittest.TestCase):
    indel = {'TYPE': 'INDEL',
             'REGION': 'INTRONIC',
             'LOCATION': 'ANY',
             'IMPACT': 5,
             "ZYGOSITY": "HETEROZYGOUS"}
    good_var = Indel(indel)
    # test 7

    def test_get_type(self):
        self.assertEqual(self.good_var.get_type(), "INDEL")
    # test 8

    def test_get_region(self):
        self.assertEqual(self.good_var.get_region(), "INTRONIC")

class InsersionMethodTesting(unittest.TestCase):
    indel_any = {'TYPE': 'INDEL',
                 'REGION': 'INTRONIC',
                 'LOCATION': 'ANY', 'IMPACT': 5,
                 "ZYGOSITY": "HETEROZYGOUS"}
    indel_any = Indel(indel_any)
    indel_spec = {'TYPE': 'INDEL',
                  'REGION': 'CODING',
                  'LOCATION': 132576250, 'IMPACT': 8,
                  "ZYGOSITY": "HETEROZYGOUS"}
    indel_spec = Indel(indel_spec)
    positive_transcript = Transcript(64805)  # SOX9
    negative_transcript = Transcript(68960)  # TOR1A

    # test 9
    def test_check_insertion_generation(self):
        self.assertEqual(len(self.indel_any.get_insertion_str(50)), 50)

    # test 10
    def test_get_insertion_any(self):
        insertion = self.indel_any.get_insertion(self.positive_transcript)
        self.assertEqual(len(insertion["alt"]), 6)

    # test 11
    def test_get_insertion_spec(self):
        insertion = self.indel_spec.get_insertion(self.negative_transcript)
        self.assertEqual(insertion["pos"], 132576250)
        self.assertEqual(len(insertion["alt"]), 9)

    # test 12
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
    # test 13

    def test_basic_deletion(self):
        indel = {'TYPE': 'INDEL',
                 'REGION': 'CODING',
                 'LOCATION': 70117532,
                 'IMPACT': -5,
                 "ZYGOSITY": "HETEROZYGOUS"}
        indel = Indel(indel)
        deletion = indel.get_deletion(self.positive_transcript)
        self.assertEqual(deletion["ref"], "ATGAAT")
        self.assertEqual(deletion["alt"], "A")
        self.assertEqual(deletion["pos"], 70117532)
    # test 14

    def test_any_location_deletion(self):
        indel = {'TYPE': 'INDEL',
                 'REGION': 'CODING',
                 'LOCATION': "ANY",
                 'IMPACT': -5,
                 "ZYGOSITY": "HETEROZYGOUS"}
        indel = Indel(indel)
        deletion = indel.get_deletion(self.positive_transcript)
        self.assertEqual(len(deletion["ref"]), 6)
        self.assertEqual(len(deletion["alt"]), 1)
    # test 15
    def test_raises_length_over_200(self):
        indel = {'TYPE': 'INDEL',
                 'REGION': 'CODING',
                 'LOCATION': "ANY",
                 'IMPACT': 570,
                 "ZYGOSITY": "HETEROZYGOUS"}
        indel = Indel(indel)
        self.assertRaises(indel.get_deletion(self.positive_transcript))
    # test 16
    def test_raises_length_greater_than_region(self):
        sox18 = Transcript(241)  # SOX18
        indel = {'TYPE': 'INDEL',
                 'REGION': 'INTRONIC',
                 'LOCATION': "ANY",
                 'IMPACT': 199,
                 "ZYGOSITY": "HETEROZYGOUS"}
        indel = Indel(indel)
        self.assertRaises(indel.get_deletion(sox18))
    # test 17
    def test_raises_length_greater_than_location(self):
        sox18 = Transcript(241)  # SOX18
        indel = {'TYPE': 'INDEL',
                 'REGION': 'INTRONIC',
                 'LOCATION': 62680315,
                 'IMPACT': 50,
                 "ZYGOSITY": "HETEROZYGOUS"}
        indel = Indel(indel)
        self.assertRaises(indel.get_deletion(sox18))


if __name__ == '__main__':
    unittest.main()
