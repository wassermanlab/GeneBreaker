from MenDelSIM.src.variant import Variant
from MenDelSIM.src.transcript import Transcript
import random


class Indel(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template: dict, transcript: Transcript):
        Variant.__init__(self, var_template, transcript)
        self.indel_amount = self.impact["INDEL_AMOUNT"]
        self.location = self.impact["LOCATION"]

        # checks
        self.check_indel()

    def check_amount(self):
        """checks amount of indel is integer between -200 and 200"""
        if type(self.indel_amount) is not int:
            raise ValueError("""For indel type the program expects an int
            representing the amount to delete or insert, 
            did not get that input.""")
        if self.indel_amount > 200 or self.indel_amount < -200:
            raise ValueError("Indels must be less than 200 in length")
        if self.indel_amount == 0:
            raise ValueError("Indel length must be not equal to 0")

    def get_region_range(self):
        """get ther region requested for that variant using the attached transcript"""
        regions = self.transcript.get_requested_region(self.region)
        region_range = []
        if (self.indel_amount > 0):  # if positive just get normal
            for region in regions:
                region_range = region_range + list(range(region[0], region[1]))
        else:
            if self.region is "CODING":
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
        if type(self.location) is int:
            self.location = self.location - 1
            if (self.indel_amount>0):
                self.check_location(self.location)
            else: 
                end = self.location + 1 - self.indel_amount
                self.check_location(self.location, end)
            return
        if self.location is not "ANY":
            raise ValueError("Location must be ANY or int. ")   

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
        if self.location == "ANY":  # pick any position within the ranges
            pos = random.choice(region_range)  # determin what this is
        else:
            pos = self.location
        return {"pos": pos,
                "ref": self.get_seq(self.transcript.get_chr(), pos, pos+1-self.indel_amount, self.transcript.get_genome()),
                "alt": self.get_seq(self.transcript.get_chr(), pos, pos+1, self.transcript.get_genome())}

    def get_insertion(self) -> dict:
        """returns (ref, alt) tuple of insersion"""
        # get requested region
        region_range = self.get_region_range()
        if self.location == "ANY":  # pick any position within the ranges
            pos = random.choice(region_range)
        else:
            pos = self.location
        return {"pos": pos,
                "ref": self.get_seq(self.transcript.get_chr(), pos, pos+1, self.transcript.get_genome()),
                "alt": self.get_seq(self.transcript.get_chr(), pos, pos+1, self.transcript.get_genome()) + self.get_insertion_str(self.indel_amount)}

    def get_vcf_row(self) -> str:
        chrom = str(self.transcript.get_chr())
        if self.indel_amount > 0:  # insersion
            var_dict = self.get_insertion()
        if self.indel_amount < 0:  # deletion
            var_dict = self.get_deletion()
        pos = str(var_dict["pos"] + 1)  # add 1 to make it on based
        ref = str(var_dict["ref"])
        alt = str(var_dict["alt"])
        ID = "_".join(["indel", pos, str(self.indel_amount)])
        if self.zygosity == "HOMOZYGOUS":
            zygosity = "1/1"
        if self.zygosity == "HEMIZYGOUS":
            zygosity = "1"
        if self.zygosity == "HETEROZYGOUS":
            zygosity = "0/1"
        return "\t".join([chrom, pos, ID, ref, alt, ".", ".", ".", "GT", zygosity])
