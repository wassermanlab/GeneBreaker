# python -m unittest tests.test_gene
import unittest
from simulator.src.variant import Variant


class VariantCreationTests(unittest.TestCase):

    def test_bad_formatting(self):
        variant = {
            "TYPE": "INDEL",
            "REGION": "INTRONIC",
            "CHECK": -8332,
            "ZYGOSITY": "HETEROZYGOUS"}
        
        with self.assertRaises(KeyError):
            Variant(variant)

    def test_region_fail(self):
        variant = {
            "TYPE": "INDEL",
            "REGION": "TEST",
            "IMPACT": {"INDEL_AMOUNT":-400 , "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        
        with self.assertRaises(ValueError):
            Variant(variant)
    
    def test_region_fail2(self):
        variant = {
            "TYPE": "INDEL",
            "REGION": "chr2:126-245a",
            "IMPACT": {"INDEL_AMOUNT":-400 , "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Variant(variant)
    
    def test_region_pass(self):
        variant = {
            "TYPE": "INDEL",
            "REGION": "chr2:126-245",
            "IMPACT": {"INDEL_AMOUNT":-400 , "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        good_var = Variant(variant)

    def test_type_fail(self):
        variant = {
            "TYPE": "TEST",
            "REGION": "INTRONIC",
            "IMPACT": {"INDEL_AMOUNT": -400, "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Variant(variant)


class VariantMethodTests(unittest.TestCase):
    good_var = {
        "TYPE": "INDEL",
        "REGION": "INTRONIC",
        "IMPACT": {"INDEL_AMOUNT": -100, "LOCATION": "ANY"},
        "ZYGOSITY": "HETEROZYGOUS"}
    good_var = Variant(good_var)

    def test_get_type(self):
        self.assertEqual(self.good_var.get_type(), "INDEL")

    def test_get_region(self):
        self.assertEqual(self.good_var.get_region(), "INTRONIC")


if __name__ == '__main__':
    unittest.main()
