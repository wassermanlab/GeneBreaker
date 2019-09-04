# python -m unittest tests.test_gene 
import unittest
from MenDelSIM.src.variants import Variants

class VariantCreationTests(unittest.TestCase):
    def test_writing_file_denovo(self):
        config = "MenDelSIM/tests/testing_data/denovo.json"
        Variants(config).save_vcf_output("MenDelSIM/tests/testing_data/denovo.output.vcf")

    def test_writing_file_biparental2(self):
        config = "MenDelSIM/tests/testing_data/biparental_2.json"
        Variants(config).save_vcf_output("MenDelSIM/tests/testing_data/biparental2.output.vcf")

if __name__ == '__main__':
    unittest.main()