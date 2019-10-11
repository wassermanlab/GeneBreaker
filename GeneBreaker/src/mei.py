from GeneBreaker.src.variant import Variant
from GeneBreaker.src.transcript import Transcript 
import random
from Bio import SeqIO
import os
from zipfile import ZipFile

class MEI(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template: dict, transcript: Transcript):
        Variant.__init__(self, var_template, transcript)
        self.element = self.impact["ELEMENT"]
        self.start = self.impact["START"]
        self.check_mei()
        # make mei
        region_range = self.get_region_range()
        if self.start == "ANY": # pick any position within the ranges
            self.pos = random.choice(region_range)
        else:
            self.pos = self.start
        self.id = "_".join(["mei", str(self.pos+1), self.element])
        self.ref = self.get_seq(self.transcript.chrom, self.pos, self.pos+1, self.transcript.genome)
        self.alt = "<INS:MEI:"+self.element+">"
        self.zygosity = var_template["ZYGOSITY"]

    def check_element(self):
        """checks element validity"""
        if self.element not in ["ALU", "LINE", "SVA"]:
            raise ValueError("""Only elements in the set of: "ALU", "LINE", "SVA" are currently supported.""")

    def check_mei(self):
        """checks all mei features"""
        if self.type != "MEI":
            raise ValueError("Must be MEI type")
        self.check_element()
        if type(self.start) == int:
            self.start = self.start - 1
            self.check_location(self.start)
        elif self.start != "ANY":
            raise ValueError("locations must be ANY or and int")

    def get_vcf_row(self) -> dict:
        element_lengths = {"ALU": 281,"LINE": 6019, "SVA": 1316}
        chrom = self.chrom
        pos = str(self.pos + 1)  # add 1 to make 1 based
        ref = str(self.ref)
        alt = str(self.alt)
        ID = str(self.id)
        if self.zygosity == "HOMOZYGOUS":
            zygosity = "1/1"
        if self.zygosity == "HEMIZYGOUS":
            zygosity = "1/1"
        if self.zygosity == "HETEROZYGOUS":
            zygosity = "0/1"
        svend = self.pos + 1 + element_lengths[self.element]
        svlen = element_lengths[self.element]
        return {
            "chrom": chrom,
            "pos":  pos,
            "id": ID,
            "ref": ref,
            "alt": alt,
            "qual": ".",
            "filter": ".",
            "info": "SVTYPE=INS;END=" + str(svend) + ";SVLEN=" + str(svlen) + ";",
            "format": "GT",
            "proband": zygosity}
