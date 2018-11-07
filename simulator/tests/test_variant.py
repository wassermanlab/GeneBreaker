# python -m unittest tests.test_gene 
import unittest
from simulator.variant import Variant

class VariantCreationTests(unittest.TestCase):

    def test_bad_formatting(self):
        variant = {
        "TYPE": "INDEL",
        "REGION": "INTRONIC",
        "CHECK": {"TYPE_IMPACT": -8332, "LOCATION": "ANY"}}
        self.assertRaises(Variant(variant))
    
    def test_region_fail(self):
        variant = {
        "TYPE": "INDEL",
        "REGION": "TEST",
        "IMPACT":{"TYPE_IMPACT": -8332, "LOCATION": "ANY"}}
        self.assertRaises(Variant(variant))

    def test_type_fail(self):
        variant = {
        "TYPE": "TEST",
        "REGION": "INTRONIC",
        "IMPACT": {"TYPE_IMPACT": -8332, "LOCATION": "ANY"}}
        self.assertRaises(Variant(variant))


class VariantMethodTests(unittest.TestCase):
    good_var = {
        "TYPE": "INDEL",
        "REGION": "INTRONIC",
        "IMPACT": {"TYPE_IMPACT": -8332, "LOCATION": "ANY"}}
    good_var = Variant(good_var)

    def test_get_type(self):
        self.assertEqual(self.good_var.get_type(), "INDEL") 

    def test_get_region(self):
        self.assertEqual(self.good_var.get_region(), "INTRONIC")

    def test_get_impact(self):
        self.assertEqual(self.good_var.get_impact(), dict({'LOCATION': 'ANY', 'TYPE_IMPACT': -8332}))

if __name__ == '__main__':
    unittest.main()