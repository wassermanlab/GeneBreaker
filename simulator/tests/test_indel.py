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
        "LOCATION": "ANY"}
        self.assertRaises(Indel(indel))
    # test 2
    def test_region_value(self):
        indel = {"TYPE": "INDEL",
        "REGION": "TEST",
        "IMPACT": 5, 
        "LOCATION": "ANY"}
        self.assertRaises(Indel(indel))
    # test 3
    def test_type_value(self):
        indel = {"TYPE": "TEST",
        "REGION": "INTRONIC",
        "IMPACT": 5,
        "LOCATION": "ANY"}
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
        "LOCATION": 0}
        self.assertRaises(Indel(indel))
    # test 6
    def test_location_str(self):
        indel = {"TYPE": "INDEL",
        "REGION": "INTRONIC",
        "IMPACT": 5, 
        "LOCATION": "four"}
        self.assertRaises(Indel(indel))    


class InheritanceVariantMethodTests(unittest.TestCase):
    indel = {'TYPE': 'INDEL',
            'REGION': 'INTRONIC',
            'LOCATION': 'ANY', 
            'IMPACT': 5}
    good_var = Indel(indel)
    # test 7
    def test_get_type(self):
        self.assertEqual(self.good_var.get_type(), "INDEL") 
    # test 8
    def test_get_region(self):
        self.assertEqual(self.good_var.get_region(), "INTRONIC")

# ========================sox9========================
# utr 5
# (70117160L, 70117532L)
# coding exons
# [(70117532L, 70117963L), (70118859L, 70119113L), (70119683L, 70120528L)]
# introns
# [(70117963L, 70118859L), (70119113L, 70119683L)]
# utr 3
# (70120528L, 70122560L)
# .[(132575220L, 132576250L), (132586364L, 132586441L)]
class InsersionMethodTesting(unittest.TestCase):
    indel_any = {'TYPE': 'INDEL',
            'REGION': 'INTRONIC',
            'LOCATION': 'ANY', 'IMPACT': 5}
    indel_any = Indel(indel_any)
    indel_spec = {'TYPE': 'INDEL',
            'REGION': 'CODING',
            'LOCATION': 132576250, 'IMPACT': 8}
    indel_spec = Indel(indel_spec)
    positive_transcript = Transcript("SOX9", 0)
    negative_transcript = Transcript("TOR1A", 0)

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
        self.assertEquals(len(insertion[4]),9)
        print insertion

class DeletionMethodTestCase(unittest.TestCase):
    positive_transcript = Transcript("SOX9", 0)
    # test 13
    def test_basic_deletion(self):
        indel = {'TYPE': 'INDEL',
            'REGION': 'CODING',
            'LOCATION': 70117532, 
            'IMPACT': -5}
        indel = Indel(indel)
        deletion = indel.get_deletion(self.positive_transcript)
        self.assertEqual(deletion["ref"], "ATGAA")
        self.assertEqual(deletion["alt"], "")
        self.assertEqual(deletion["pos"], 70117532)
    # test 14
    def test_any_location_deletion(self):
        indel = {'TYPE': 'INDEL',
            'REGION': 'CODING',
            'LOCATION': "ANY", 
            'IMPACT': -5}
        indel = Indel(indel)
        deletion = indel.get_deletion(self.positive_transcript)
        self.assertEqual(len(deletion["ref"]), 5)
        self.assertEqual(deletion["alt"], "")
    # test 15
    def test_raises_length_over_200(self):
        indel = {'TYPE': 'INDEL',
            'REGION': 'CODING',
            'LOCATION': "ANY", 
            'IMPACT': 570}
        indel = Indel(indel)
        self.assertRaises(indel.get_deletion(self.positive_transcript))
    # test 16
    def test_raises_length_greater_than_region(self):
        sox18 = Transcript("SOX18", 0)
        indel = {'TYPE': 'INDEL',
            'REGION': 'INTRONIC',
            'LOCATION': "ANY", 
            'IMPACT': 199}
        indel = Indel(indel)
        self.assertRaises(indel.get_deletion(sox18))
    # test 17
    def test_raises_length_greater_than_location(self):
        sox18 = Transcript("SOX18", 0)
        indel = {'TYPE': 'INDEL',
            'REGION': 'INTRONIC',
            'LOCATION': 62680315, 
            'IMPACT': 50}
        indel = Indel(indel)
        self.assertRaises(indel.get_deletion(sox18))


if __name__ == '__main__':
    unittest.main()