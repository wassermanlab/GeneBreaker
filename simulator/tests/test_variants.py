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

    def test_writing_file_paternal(self):
        config = "simulator/tests/testing_data/paternal.json"
        Variants(config).save_vcf_output("simulator/tests/testing_data/paternal.output.vcf")

    def test_writing_file_maternal(self):
        config = "simulator/tests/testing_data/maternal.json"
        Variants(config).save_vcf_output("simulator/tests/testing_data/maternal.output.vcf")

    def test_writing_file_denovo(self):
        config = "simulator/tests/testing_data/denovo.json"
        Variants(config).save_vcf_output("simulator/tests/testing_data/denovo.output.vcf")

    def test_writing_file_biparental1(self):
        config = "simulator/tests/testing_data/biparental_1.json"
        Variants(config).save_vcf_output("simulator/tests/testing_data/biparental1.output.vcf")
    
    def test_writing_file_biparental2(self):
        config = "simulator/tests/testing_data/biparental_2.json"
        Variants(config).save_vcf_output("simulator/tests/testing_data/biparental2.output.vcf")

if __name__ == '__main__':
    unittest.main()