# python -m unittest tests.test_gene 
import unittest
from simulator.short_tandem_repeat import ShortTandemRepeat

class STRBasicTests(unittest.TestCase):
    # test 1
    def test_wrong_keys(self):
        STR = {
        "IMPACT": {
            "CHROM": "chr2", 
            "END": 191745646, 
            "START": 191745599, 
            "STR": 5
        }, 
        "LOCATION": "NONE", 
        "REGION": "UTR", 
        "TYPE": "STR"}
        STR = ShortTandemRepeat(STR)
        print STR.get_str()


if __name__ == '__main__':
    unittest.main()