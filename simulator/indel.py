from simulator.variant import Variant
from simulator.transcript import Transcript 
import random

class Indel(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template):
        try:
            Variant.__init__(self, var_template)
            if type(self.impact) is not int:
                raise Exception("""For indel type the program expects an int
                representing the amount to delete or insert, 
                did not get that input.""")
            if self.impact > 200 or self.impact < -200:
                raise Exception("Indels must be less than 200")
            if self.impact == 0:
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
    
    def get_deletion(self, transcript):
        """returns (pos ,ref, alt) tuple of deletion"""
        # get requested region
        regions = transcript.get_requested_region(self.region)
        if len(regions) == 0:
            raise Exception("region requested does not exists")
        # get ranges
        region_range = []
        for region in regions: # range must be cut so that there is no overlap
            region_range = region_range + range(region[0], region[1]+self.impact)
        if len(region_range) == 0:
            raise Exception("""regions selected are too small to accomidate a 
            deletion of this size, try reducing the size of the deletion""")
        if self.location == "ANY": # pick any position within the ranges
            pos = random.choice(region_range)  # determin what this is 
            shifted_pos = pos - transcript.get_start()
            return {"pos": pos,
                    "ref": transcript.get_seq()[shifted_pos:shifted_pos-self.impact],
                    "alt": ""}
        else:
            # check that deletion doesn't go over the amount
            if self.location not in region_range:
                raise Exception("""location selected is to0 small to accomidate
                a deletion of this size, try reducing the size of the 
                deletion""")
            pos = self.location
            shifted_pos = pos - transcript.get_start()
            return {"pos": pos,
                    "ref": transcript.get_seq()[shifted_pos:shifted_pos-self.impact],
                    "alt": ""}

    def get_insertion(self, transcript):
        """returns (ref, alt) tuple of insersion"""
        # get requested region
        regions = transcript.get_requested_region(self.region)
        if len(regions) == 0:
            raise Exception("region requested does not exists")
        if self.location == "ANY": # pick any position within the ranges
            # get ranges
            region_range = []
            for region in regions:
                region_range = region_range + range(region[0], region[1])
            pos = random.choice(region_range)
            shifted_pos = pos - transcript.get_start()
            return {"pos": pos,
                    "ref": transcript.get_seq()[shifted_pos],
                    "alt": transcript.get_seq()[shifted_pos] + self.get_insertion_str(self.impact)}
        else:
            pos = self.location
            shifted_pos = pos - transcript.get_start()
            return {"pos": pos,
                    "ref": transcript.get_seq()[shifted_pos],
                    "alt": transcript.get_seq()[shifted_pos] + self.get_insertion_str(self.impact)}

    def get_vcf_row(self, transcript):
        chrom = str(transcript.get_chr())
        if self.impact > 0:  # insersion
            var_dict = self.get_insertion(transcript)
        if self.impact < 0: # deletion
            var_dict = self.get_deletion(transcript)
        pos = str(var_dict["pos"] + 1) # add 1 to make itone based
        ref = str(var_dict["ref"])
        alt = str(var_dict["alt"])
        ID = "_".join(["indel", pos, str(self.impact)])
        return "\t".join([chrom, pos, ID, ref, alt])
