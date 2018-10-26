import unittest
from simulator.gene import Gene

class GeneSetup(unittest.TestCase):
    def setUp(self):
        self.simple_gene = Gene("testing_data/genes/basic_gene.json")
    
    def tearDown(self):
        self.simple_gene = None


# class GeneMethodTests(GeneSetup):
#     def get_exons(self):
#         return False

#     def get_introns(self):
#         return False

#     def get_upstream_utr(self):
#         return False

#     def get_downstream_utr(self):
#         return False

#     def get_promoter(self):
#         return False

#     def get_enhancer(self):
        

if __name__ == '__main__':
    unittest.main()