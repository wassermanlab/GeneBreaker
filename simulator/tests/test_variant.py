# python -m unittest tests.test_gene 
import unittest
from simulator.variant import Variant

class VariantCreationTests(unittest.TestCase):

    def test_bad_formatting(self):
        variant = {
        "TYPE": "INDEL",
        "REGION": "INTRONIC",
        "CHECK": -8332, 
        "LOCATION": "ANY"}
        self.assertRaises(Variant(variant))
    
    def test_region_fail(self):
        variant = {
        "TYPE": "INDEL",
        "REGION": "TEST",
        "IMPACT": -8332, 
        "LOCATION": "ANY"}
        self.assertRaises(Variant(variant))

    def test_type_fail(self):
        variant = {
        "TYPE": "TEST",
        "REGION": "INTRONIC",
        "IMPACT": -8332, 
        "LOCATION": "ANY"}
        self.assertRaises(Variant(variant))


class VariantMethodTests(unittest.TestCase):
    good_var = {
        "TYPE": "INDEL",
        "REGION": "INTRONIC",
        "IMPACT": -8332, 
        "LOCATION": "ANY"}
    good_var = Variant(good_var)

    def test_get_type(self):
        self.assertEqual(self.good_var.get_type(), "INDEL") 

    def test_get_region(self):
        self.assertEqual(self.good_var.get_region(), "INTRONIC")

if __name__ == '__main__':
    unittest.main()