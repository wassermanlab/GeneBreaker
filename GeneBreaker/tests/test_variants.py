# python -m unittest tests.test_gene 
import unittest
from GeneBreaker.src.variants import Variants

class VariantCreationTests(unittest.TestCase):
    def test_writing_file_denovo(self):
        config = "GeneBreaker/tests/testing_data/denovo.json"
        _vars = Variants(config).get_variant_rows()

    def test_writing_file_biparental2(self):
        config = "GeneBreaker/tests/testing_data/biparental_2.json"
        _vars = Variants(config).get_variant_rows()



if __name__ == '__main__':
    unittest.main()