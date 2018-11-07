# python -m unittest tests.test_gene 
import unittest
from simulator.variants import Variants

class VariantCreationTests(unittest.TestCase):

    def test_bad_formatting(self):
        gene = "simulator/tests/testing_data/genes/basic_gene.json"
        variants = "simulator/tests/testing_data/variants/vars1.json"
        self.assertRaises(Variants(gene, variants))
    
    def test_good_formatting(self):
        gene = "simulator/tests/testing_data/genes/basic_gene.json"
        variants = "simulator/tests/testing_data/variants/vars2.json"
        print Variants(gene, variants).variants_2_VCF()
    
    def test_writing_file(self):
        gene = "simulator/tests/testing_data/genes/basic_gene.json"
        variants = "simulator/tests/testing_data/variants/vars2.json"
        Variants(gene, variants).save_vcf_output("output.vcf")


if __name__ == '__main__':
    unittest.main()