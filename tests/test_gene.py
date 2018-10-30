# python -m unittest tests.test_gene 
import unittest
from simulator.gene import Gene

class GeneSetup(unittest.TestCase):
    def setUp(self):
        self.simple_gene = Gene("tests/testing_data/genes/basic_gene.json")
    
    def tearDown(self):
        self.simple_gene = None

class GeneCreationTests(unittest.TestCase): 
    def test_IO_exception(self):
        self.assertRaises(Gene("tests/genes/basic_gene_nonexist.json"))
    
    def test_bad_json_exception(self):
        self.assertRaises(Gene("tests/testing_data/genes/impropper_format_gene.json"))

class GeneMethodTests(GeneSetup):

    def test_get_exons(self):
        self.assertEqual(self.simple_gene.get_exons(), 
                        [[0, 96],
                        [8333, 8902]])
        self.assertEqual(self.simple_gene.get_requested_region("CODING"), 
                        [[0, 96],
                        [8333, 8902]])

    def test_get_introns(self):
        self.assertEqual(self.simple_gene.get_introns(), [[96,8333]])
        self.assertEqual(self.simple_gene.get_requested_region("INTRONIC"), [[96,8333]])


    # def test_get_utr(self):
    #     self.assertTrue(self.simple_gene.get_utr())

    # def test_get_promoter(self):
    #     self.assertTrue(self.simple_gene.get_promoter())

    # def test_get_enhancer(self):
    #     self.assertTrue(self.simple_gene.get_enhancer())


if __name__ == '__main__':
    unittest.main()