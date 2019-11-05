# python -m unittest tests.test_gene
import unittest
from GeneBreaker.src.variant import Variant
from GeneBreaker.src.transcript import Transcript
from GeneBreaker.src.api_helper import *

class VariantCreationTests(unittest.TestCase):
    XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
    transcript = Transcript(XKR8_uid, "hg38")

    def test_bad_formatting(self):
        variant = {
            "TYPE": "INDEL",
            "REGION": "INTRONIC",
            "CHECK": -8332,
            "ZYGOSITY": "HETEROZYGOUS"}
        
        with self.assertRaises(KeyError):
            Variant(variant, self.transcript)

    def test_region_fail(self):
        variant = {
            "TYPE": "INDEL",
            "REGION": "TEST",
            "IMPACT": {"INDEL_AMOUNT":-400 , "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        
        with self.assertRaises(ValueError):
            Variant(variant, self.transcript)
    
    def test_region_fail2(self):
        variant = {
            "TYPE": "INDEL",
            "REGION": "2:126-245a",
            "IMPACT": {"INDEL_AMOUNT":-400 , "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Variant(variant, self.transcript)
    
    def test_region_pass(self):
        variant = {
            "TYPE": "INDEL",
            "REGION": "2:126-245",
            "IMPACT": {"INDEL_AMOUNT":-400 , "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        good_var = Variant(variant, self.transcript)

    def test_type_fail(self):
        variant = {
            "TYPE": "TEST",
            "REGION": "INTRONIC",
            "IMPACT": {"INDEL_AMOUNT": -400, "LOCATION": "ANY"},
            "ZYGOSITY": "HETEROZYGOUS"}
        with self.assertRaises(ValueError):
            Variant(variant, self.transcript)


class VariantMethodTests(unittest.TestCase):
    XKR8_uid = get_all_transcripts("XKR8", "hg38")[0]["qualifiers"]["uid"]
    transcript = Transcript(XKR8_uid, "hg38")
    
    good_var = {
        "TYPE": "INDEL",
        "REGION": "INTRONIC",
        "IMPACT": {"INDEL_AMOUNT": -100, "LOCATION": "ANY"},
        "ZYGOSITY": "HETEROZYGOUS"}
    good_var = Variant(good_var, transcript)

    def test_get_type(self):
        self.assertEqual(self.good_var.get_type(), "INDEL")

    def test_get_region(self):
        self.assertEqual(self.good_var.get_region(), "INTRONIC")


if __name__ == '__main__':
    unittest.main()
