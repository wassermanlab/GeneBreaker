# python -m unittest tests.test_gene 
import unittest
from simulator.indel import Indel
import json

class IndelCreationTests(unittest.TestCase):
    f = open("tests/testing_data/variants/vars8.json")
    wrong_keys = json.load(f)
    wrong_keys = wrong_keys["VAR1"]
    f = open("tests/testing_data/variants/vars9.json")
    region_fail = json.load(f)
    region_fail = region_fail["VAR1"]
    f = open("tests/testing_data/variants/vars10.json")
    type_fail = json.load(f)
    type_fail = type_fail["VAR1"]
    f = open("tests/testing_data/variants/vars11.json")
    type_fail2 = json.load(f)
    type_fail2 = type_fail2["VAR1"]

    def test_bad_formatting(self):
        self.assertRaises(Indel(self.wrong_keys))
    
    def test_region_fail(self):
        self.assertRaises(Indel(self.region_fail))

    def test_type_fail(self):
        self.assertRaises(Indel(self.type_fail))

    def test_type_fail2(self):
        self.assertRaises(Indel(self.type_fail2))


class InheritanceVariantMethodTests(unittest.TestCase):
    f = open("tests/testing_data/variants/vars6.json")
    good_var = json.load(f)
    good_var = good_var["VAR1"]
    good_var = Indel(good_var)

    def test_get_type(self):
        self.assertEqual(self.good_var.get_type(), "INDEL") 

    def test_get_region(self):
        self.assertEqual(self.good_var.get_region(), "INTRONIC")

    def test_get_impact(self):
        self.assertEqual(self.good_var.get_impact(), 5)

if __name__ == '__main__':
    unittest.main()