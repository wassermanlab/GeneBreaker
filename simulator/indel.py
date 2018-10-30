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
                representing the amount to delete or insert, 
                did not get that input.""")
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
            # randomly choses ACGT for each position in insertion
            insertion = insertion + random.choice("ACGT")
        return insertion

    def get_deletion(self, gene):
        """returns (pos ,ref, alt) tuple of deletion"""
        # get requested region
        regions = gene.get_requested_region(self.region)
        if len(regions) == 0:
            raise Exception("region requested does not exists")
        # get ranges
        region_range = []
        for region in regions: # range must be cut so that there is no overlap
            region_range = region_range + \
            range(region[0], region[1]+self.type_impact)
        if len(region_range) == 0:
            raise Exception("""regions selected are to0 small to accomidate a 
            deletion of this size, try reducing the size of the deletion""")
        if self.location == "ANY": # pick any position within the ranges
            pos = random.choice(region_range)  # determin what this is 
            return {"pos": pos,
                    "ref": gene.get_seq()[pos:pos-self.type_impact],
                    "alt": ""}
        else:
            # check that deletion doesn't go over the amount
            if self.location not in region_range:
                raise Exception("""location selected is to0 small to accomidate
                a deletion of this size, try reducing the size of the 
                deletion""")
            return {"pos": self.location,
                    "ref": gene.get_seq()[self.location:self.location-self.type_impact],
                    "alt": ""}

    def get_insertion(self, gene):
        """returns (ref, alt) tuple of insersion"""
        # get requested region
        regions = gene.get_requested_region(self.region)
        if len(regions) == 0:
            raise Exception("region requested does not exists")
        if self.location == "ANY": # pick any position within the ranges
            # get ranges
            region_range = []
            for region in regions:
                region_range = region_range + range(region[0], region[1])
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
        if self.type_impact > 0:  # insersion
            var_dict = self.get_insertion(gene)
        if self.type_impact < 0: # deletion
            var_dict = self.get_deletion(gene)
        pos = str(var_dict["pos"])
        ref = str(var_dict["ref"])
        alt = str(var_dict["alt"])
        ID = "_".join(["indel", pos, str(self.type_impact)])
        return "\t".join([chrom, pos, ID, ref, alt])
