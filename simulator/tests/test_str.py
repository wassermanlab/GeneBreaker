# python -m unittest tests.test_gene 
import unittest
from simulator.short_tandem_repeat import ShortTandemRepeat
from simulator.transcript import Transcript

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
        self.assertEquals(STR.get_str_motif(), "GCA")
    
    def test_retraction_too_large(self):
        STR = {
        "IMPACT": {
            "CHROM": "chr2", 
            "END": 191745646, 
            "START": 191745599, 
            "STR": -50
        }, 
        "LOCATION": "NONE", 
        "REGION": "UTR", 
        "TYPE": "STR"}
        STR = ShortTandemRepeat(STR)
        with self.assertRaises(Exception) as cm:
            STR.get_retraction()
        err = cm.exception
        self.assertEqual(str(err), 'retraction length is larger than the total str')
    
    def test_retraction_works(self):
        STR = {
        "IMPACT": {
            "CHROM": "chr2", 
            "END": 191745646, 
            "START": 191745599, 
            "STR": -5
        }, 
        "LOCATION": "NONE", 
        "REGION": "UTR", 
        "TYPE": "STR"}
        STR = ShortTandemRepeat(STR)
        retraction = STR.get_retraction() 
        self.assertEquals(retraction['pos'], 191745599)
        self.assertEquals(retraction['ref'], "GCAGCAGCAGCAGCA")
        self.assertEquals(retraction['alt'], "")

    def test_insertion_too_large(self):
        STR = {
        "IMPACT": {
            "CHROM": "chr2", 
            "END": 191745646, 
            "START": 191745599, 
            "STR": 10000
        }, 
        "LOCATION": "NONE", 
        "REGION": "UTR", 
        "TYPE": "STR"}
        STR = ShortTandemRepeat(STR)
        with self.assertRaises(Exception) as cm:
            STR.get_expantion()
        err = cm.exception
        self.assertEqual(str(err), 'expansion length is too large') 

    def test_insertion_works(self):
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
        expansion =  STR.get_expantion()
        self.assertEquals(expansion['pos'], 191745599)
        self.assertEquals(expansion['ref'], "G")
        self.assertEquals(expansion['alt'], "GCAGCAGCAGCAGCAG")

if __name__ == '__main__':
    unittest.main()