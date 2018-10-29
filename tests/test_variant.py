# python -m unittest tests.test_gene 
import unittest
from simulator.variant import Variant
import json

class VariantCreationTests(unittest.TestCase):
    f = open("tests/testing_data/variants/vars8.json")
    wrong_keys = json.load(f)
    wrong_keys = wrong_keys["VAR1"]
    f = open("tests/testing_data/variants/vars9.json")
    region_fail = json.load(f)
    region_fail = region_fail["VAR1"]
    f = open("tests/testing_data/variants/vars10.json")
    type_fail = json.load(f)
    type_fail = type_fail["VAR1"]

    def test_bad_formatting(self):
        self.assertRaises(Variant(self.wrong_keys))
    
    def test_region_fail(self):
        self.assertRaises(Variant(self.region_fail))

    def test_type_fail(self):
        self.assertRaises(Variant(self.type_fail))


class VariantMethodTests(unittest.TestCase):
    f = open("tests/testing_data/variants/vars7.json")
    good_var = json.load(f)
    good_var = good_var["VAR1"]
    good_var = Variant(good_var)

    def test_get_type(self):
        self.assertEqual(self.good_var.get_type(), "INDEL") 

    def test_get_region(self):
        self.assertEqual(self.good_var.get_region(), "INTRONIC")

    def test_get_impact(self):
        self.assertEqual(self.good_var.get_impact(), -8332)

if __name__ == '__main__':
    unittest.main()