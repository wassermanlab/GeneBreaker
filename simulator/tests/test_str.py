# python -m unittest tests.test_gene 
import unittest
from simulator.short_tandem_repeat import ShortTandemRepeat

class STRBasicTests(unittest.TestCase):
    # test 1
    def test_correct_motif(self):
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
        print STR.get_str_motif()
    
    # def test_retraction_too_large(self):
    #     positive_transcript = Transcript("GLS", 0)
    #     STR = {
    #     "IMPACT": {
    #         "CHROM": "chr2", 
    #         "END": 191745646, 
    #         "START": 191745599, 
    #         "STR": -50
    #     }, 
    #     "LOCATION": "NONE", 
    #     "REGION": "UTR", 
    #     "TYPE": "STR"}
    #     STR = ShortTandemRepeat(STR)
    #     self.assertRaises(STR.get_retraction())
    
    # def test_retraction_works(self):
    #     positive_transcript = Transcript("GLS", 0)
    #     STR = {
    #     "IMPACT": {
    #         "CHROM": "chr2", 
    #         "END": 191745646, 
    #         "START": 191745599, 
    #         "STR": -5
    #     }, 
    #     "LOCATION": "NONE", 
    #     "REGION": "UTR", 
    #     "TYPE": "STR"}
    #     STR = ShortTandemRepeat(STR)
    #     print STR.get_retraction()

    # def test_insertion_too_large(self):
    #     positive_transcript = Transcript("GLS", 0)
    #     STR = {
    #     "IMPACT": {
    #         "CHROM": "chr2", 
    #         "END": 191745646, 
    #         "START": 191745599, 
    #         "STR": 10000
    #     }, 
    #     "LOCATION": "NONE", 
    #     "REGION": "UTR", 
    #     "TYPE": "STR"}
    #     STR = ShortTandemRepeat(STR)
    #     self.assertRaises(STR.get_expantion())

    # def test_insertion_works(self):
    #     positive_transcript = Transcript("GLS", 0)
    #     STR = {
    #     "IMPACT": {
    #         "CHROM": "chr2", 
    #         "END": 191745646, 
    #         "START": 191745599, 
    #         "STR": 5
    #     }, 
    #     "LOCATION": "NONE", 
    #     "REGION": "UTR", 
    #     "TYPE": "STR"}
    #     STR = ShortTandemRepeat(STR)
    #     print STR.get_expantion()

    # def test_full(self):
    #     None
if __name__ == '__main__':
    unittest.main()