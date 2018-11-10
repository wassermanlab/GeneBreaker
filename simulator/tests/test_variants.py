# python -m unittest tests.test_gene 
import unittest
from simulator.variants import Variants

class VariantCreationTests(unittest.TestCase):
    
    def test_good_formatting(self):
        config = "simulator/tests/testing_data/test.json"
        print Variants(config).variants_2_VCF()
    
    def test_writing_file(self):
        config = "simulator/tests/testing_data/test.json"
        Variants(config).save_vcf_output("simulator/tests/testing_data/output.vcf")

if __name__ == '__main__':
    unittest.main()