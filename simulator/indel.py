from simulator.variant import Variant
from simulator.gene import Gene
import random 

class Indel(Variant): 
    # assume var_template is of type dict already
    def __init__(self, var_template):
        try: 
            Variant.__init__(self, var_template)
            if type(self.type_impact) is not int:
                raise Exception("""For indel type the program expects an int 
                representing the amount to delete or insert, did not get that input.""")
            if self.type_impact == 0:
                raise Exception("Indel length must be not equal to 0")
            if self.type != "INDEL":
                raise Exception("Must be indel type")
        except:
            print('check that type is indel and that impact is an int')
    
    def get_insertion_str(self, size):
        """returns a random string of ACGT according to size"""
        insertion = ""
        for i in range(size):
            insertion = insertion + random.choice("ACGT") # randomly choses ACGT for each position in insertion 
        return insertion

    def get_deletion(self, gene):
        """returns (pos ,ref, alt) tuple of deletion"""
        None
        # # get requested region
        # regions = gene.get_requested_region(self.region)
        # # get ranges
        # region_range = []
        # for region in regions:
        #     region_range = region_range + range(region[1], region[2])     
        # ## checks

    def get_insertion(self, gene): 
        """returns (ref, alt) tuple of insersion"""
                # get requested region
        regions = gene.get_requested_region(self.region)
        # get ranges
        region_range = []
        for region in regions:
            region_range = region_range + range(region[0], region[1]) 
        if self.location == "ANY": # pick any position within the ranges 
            pos = random.choice(region_range)
            return {"pos": pos, 
                    "ref": gene.get_seq()[pos], 
                    "alt": gene.get_seq()[pos] + self.get_insertion_str(self.type_impact)}
        else: 
            return {"pos": self.location, 
                    "ref": gene.get_seq()[self.location], 
                    "alt": gene.get_seq()[self.location] + self.get_insertion_str(self.type_impact)}

    def get_vcf_row(self, gene):
        chrom = str(gene.get_chr()) 
        if self.type_impact > 0: # insersion
            var_dict = self.get_insertion(gene)
        if self.type_impact < 0: # deletion 
            var_dict = self.get_deletion(gene)
        pos = var_tuple["pos"]
        ref = var_tuple["ref"]
        alt = var_tuple["alt"]
        ID = "_".join(["indel", pos, self.type_impact])
        return "\t".join([chr, pos, ID, ref, alt])
