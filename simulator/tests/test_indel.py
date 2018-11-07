# python -m unittest tests.test_gene 
import unittest
from simulator.indel import Indel
from simulator.gene import Gene

class IndelCreationTests(unittest.TestCase):
    # test 1
    def test_wrong_keys(self):
        indel = {
        "TYPE": "INDEL",
        "REGION": "INTRONIC",
        "CHECK": {"TYPE_IMPACT": -8332, "LOCATION": "ANY"}}
        self.assertRaises(Indel(indel))
    # test 2
    def test_region_value(self):
        indel = {"TYPE": "INDEL",
        "REGION": "TEST",
        "IMPACT":{"TYPE_IMPACT": -8332, "LOCATION": "ANY"}}
        self.assertRaises(Indel(indel))
    # test 3
    def test_type_value(self):
        indel = {"TYPE": "TEST",
        "REGION": "INTRONIC",
        "IMPACT": {"TYPE_IMPACT": -8332, "LOCATION": "ANY"}}
        self.assertRaises(Indel(indel))
    # test 4
    def test_type_value_indel(self):
        indel = {"TYPE": "SNV",
        "REGION": "INTRONIC",
        "IMPACT": {"TYPE_IMPACT": -8332, "LOCATION": "ANY"}}
        self.assertRaises(Indel(indel))
    # test 5
    def test_location_0(self):
        indel = {"TYPE": "INDEL",
        "REGION": "INTRONIC",
        "IMPACT": {"TYPE_IMPACT": -8332, "LOCATION": 0}}
        self.assertRaises(Indel(indel))
    # test 6
    def test_location_str(self):
        indel = {"TYPE": "INDEL",
        "REGION": "INTRONIC",
        "IMPACT": {"TYPE_IMPACT": -8332, "LOCATION": "four"}}
        self.assertRaises(Indel(indel))    


class InheritanceVariantMethodTests(unittest.TestCase):
    indel = {'TYPE': 'INDEL',
            'REGION': 'INTRONIC',
            'IMPACT': {'LOCATION': 'ANY', 'TYPE_IMPACT': 5}}
    good_var = Indel(indel)
    # test 7
    def test_get_type(self):
        self.assertEqual(self.good_var.get_type(), "INDEL") 
    # test 8
    def test_get_region(self):
        self.assertEqual(self.good_var.get_region(), "INTRONIC")
    # test 9
    def test_get_impact(self):
        self.assertEqual(self.good_var.get_impact(), dict({'LOCATION': 'ANY', 'TYPE_IMPACT': 5}))


class InsersionMethodTesting(unittest.TestCase):
    indel_any = {'TYPE': 'INDEL',
            'REGION': 'INTRONIC',
            'IMPACT': {'LOCATION': 'ANY', 'TYPE_IMPACT': 5}}
    indel_any = Indel(indel_any)
    simple_gene = Gene("simulator/tests/testing_data/genes/basic_gene.json")
    indel_spec = {'TYPE': 'INDEL',
            'REGION': 'CODING',
            'IMPACT': {'LOCATION': 5, 'TYPE_IMPACT': 8}}
    indel_spec = Indel(indel_spec)
    
    # test 10
    def test_check_insertion_generation(self):
        self.assertEqual(len(self.indel_any.get_insertion_str(50)), 50)
    
    # test 11
    def test_get_insertion_any(self):
        insertion = self.indel_any.get_insertion(self.simple_gene)
        self.assertEqual(len(insertion["alt"]), 6)
    
    # test 12
    def test_get_insertion_spec(self):
        insertion = self.indel_spec.get_insertion(self.simple_gene)
        self.assertEqual(insertion["pos"], 5)
    
    # test 13
    def test_get_insertion_row(self):
        insertion = self.indel_spec.get_vcf_row(self.simple_gene)

class DeletionMethodTestCase(unittest.TestCase):
    simple_gene = Gene("simulator/tests/testing_data/genes/basic_gene.json")
    # test 14
    def test_basic_deletion(self):
        indel = {'TYPE': 'INDEL',
            'REGION': 'CODING',
            'IMPACT': {'LOCATION': 5, 'TYPE_IMPACT': -5}}
        indel = Indel(indel)
        deletion = indel.get_deletion(self.simple_gene)
        self.assertEqual(deletion["ref"], "TCAGG")
        self.assertEqual(deletion["alt"], "")
        self.assertEqual(deletion["pos"], 5)
    # test 15
    def test_any_location_deletion(self):
        indel = {'TYPE': 'INDEL',
            'REGION': 'CODING',
            'IMPACT': {'LOCATION': "ANY", 'TYPE_IMPACT': -5}}
        indel = Indel(indel)
        deletion = indel.get_deletion(self.simple_gene)
        self.assertEqual(len(deletion["ref"]), 5)
        self.assertEqual(deletion["alt"], "")
    # test 16
    def test_raises_length_any(self):
        indel = {'TYPE': 'INDEL',
            'REGION': 'CODING',
            'IMPACT': {'LOCATION': "ANY", 'TYPE_IMPACT': 570}}
        indel = Indel(indel)
        self.assertRaises(indel.get_deletion(self.simple_gene))
    # test 17
    def test_raises_length_specific(self):
        indel = {'TYPE': 'INDEL',
            'REGION': 'CODING',
            'IMPACT': {'LOCATION': 5, 'TYPE_IMPACT': 96}}
        indel = Indel(indel)
        self.assertRaises(indel.get_deletion(self.simple_gene))


if __name__ == '__main__':
    unittest.main()