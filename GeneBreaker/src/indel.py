from GeneBreaker.src.variant import Variant
from GeneBreaker.src.transcript import Transcript
import random

## TODO: change indel_amount to length
class Indel(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template: dict, transcript: Transcript):
        try:
            Variant.__init__(self, var_template, transcript)
            self.indel_amount = self.impact["INDEL_AMOUNT"]
            self.pos = self.impact["START"]
            self.check_indel()
            if self.indel_amount > 0:  # get insersion
                var_dict = self.get_insertion()
            if self.indel_amount < 0:  # get deletion
                var_dict = self.get_deletion()
            self.id = "_".join(["indel", str(self.pos), str(self.indel_amount)])
            self.ref = str(var_dict["ref"])
            self.alt = str(var_dict["alt"])
        except Exception as e:
            raise(e)

    def check_amount(self):
        """checks amount of indel is integer between -200 and 200"""
        if type(self.indel_amount) != int:
            raise ValueError("""For indel type the program expects an int
            representing the amount to delete or insert, 
            did not get that input.""")
        if self.indel_amount > 200 or self.indel_amount < -200:
            raise ValueError("Indels must be less than 200 in length")
        if self.indel_amount == 0:
            raise ValueError("Indel length must be not equal to 0")

    def get_region_range(self):
        """get the region requested for that variant using the attached transcript"""
        regions = self.transcript.get_requested_region(self.region)
        region_range = []
        if (self.indel_amount > 0):  # if positive just get normal
            for region in regions:
                region_range = region_range + list(range(region[0], region[1]))
        else:
            if self.region == "CODING":
                for region in regions:
                    region_range = region_range + \
                        list(range(region[0]+self.indel_amount+1, region[1]))
            else:
                for region in regions:
                    region_range = region_range + \
                        list(range(region[0], region[1]+self.indel_amount))
        if len(region_range) == 0:
            raise ValueError("""regions selected are too small to accommodate a 
            deletion of this size, try reducing the size of the deletion""")
        return region_range

    def check_indel(self):
        """checks all indel features"""
        if self.type != "INDEL":
            raise ValueError("Must be INDEL type")
        self.check_amount()
        if type(self.pos) == int:
            self.pos = self.pos - 1
            if (self.indel_amount>0):
                self.check_location(self.pos)
            else: 
                end = self.pos + 1 - self.indel_amount
                self.check_location(self.pos, end)
            return
        if self.pos != "ANY":
            raise ValueError("Location must be ANY or int.")   

    def get_insertion_str(self, size: int) -> str:
        """returns a random string of ACGT according to size"""
        insertion = ""
        for i in list(range(size)):
            # randomly choses ACGT for each position in insertion
            insertion = insertion + random.choice("ACGT")
        return insertion

    def get_deletion(self) -> dict:
        """returns (pos ,ref, alt) tuple of deletion"""
        # get ranges
        region_range = self.get_region_range()
        if self.pos == "ANY":  # pick any position within the ranges
            pos = random.choice(region_range)  # determin what this is
        else:
            pos = self.pos
        return {"pos": pos,
                "ref": self.get_seq(self.transcript.get_chr(), pos, pos+1-self.indel_amount, self.transcript.get_genome()),
                "alt": self.get_seq(self.transcript.get_chr(), pos, pos+1, self.transcript.get_genome())}

    def get_insertion(self) -> dict:
        """returns (ref, alt) tuple of insersion"""
        # get requested region
        region_range = self.get_region_range()
        if self.pos == "ANY":  # pick any position within the ranges
            pos = random.choice(region_range)
        else:
            pos = self.pos
        return {"pos": pos,
                "ref": self.get_seq(self.transcript.get_chr(), pos, pos+1, self.transcript.get_genome()),
                "alt": self.get_seq(self.transcript.get_chr(), pos, pos+1, self.transcript.get_genome()) + self.get_insertion_str(self.indel_amount)}

  