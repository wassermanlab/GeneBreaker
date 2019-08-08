from MenDelSIM.src.variant import Variant
from MenDelSIM.src.transcript import Transcript
import random
from MenDelSIM.src.api_helper import *
from lxml import etree


class ShortTandemRepeat(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template, transcript: Transcript):
        Variant.__init__(self, var_template, transcript)
        self.chrom = self.impact["CHROM"]
        self.start = self.impact["START"]
        self.end = self.impact["END"]
        self.length = self.impact["STR"]
        self.check_str()
    # TODO: check location for retractions UTR and INTRONS
    def check_str(self):
        if self.type != "STR":
            raise ValueError("Must be STR type")
        if self.length == 0:
            raise Exception("STR length must be not equal to 0")
        self.check_location(self.start)
        self.start = self.start - 1

    def get_str_motif(self) -> str:
        """return the motif of the specific str"""
        STR = get_str(self.start+1, self.end, self.chrom,
                      self.transcript.get_genome())
        seq = STR['qualifiers']['motif']
        return seq.upper()

    def get_anchor_position(self) -> str:
        """Retrieve a DNA sequence from UCSC.
        Note: UCSC assumes 1 based indexing so we add a 1"""
        seq = self.get_seq(self.chrom, (self.start-1), self.start, self.transcript.get_genome())
        return seq.upper()

    def get_retraction(self) -> dict:
        """returns (pos ,ref, alt) tuple of retraction"""
        # get str
        STR = self.get_str_motif()
        size = len(STR)
        total_repeat_length = self.end - self.start
        # check that requested size is not over
        positive_len = (self.length*-1)
        if total_repeat_length < positive_len*size:
            raise Exception("retraction length is larger than the total str")
        alt = self.get_anchor_position()
        ref = alt + STR*positive_len
        return {"pos": self.start-1,
                "ref": ref,
                "alt": alt}

    def get_expantion(self) -> dict:
        """returns (pos ,ref, alt) tuple of retraction"""
        STR = self.get_str_motif()
        size = len(STR)
        if size*self.length > 20000:
            raise Exception("expansion length is too large")
        ref = STR
        alt = STR*self.length + STR
        return {"pos": self.start,
                "ref": ref,
                "alt": alt}

    def get_vcf_row(self) -> str:
        # get regions
        if self.length > 0:  # insersion
            var_dict = self.get_expantion()
        if self.length < 0:  # deletion
            var_dict = self.get_retraction()
        chrom = transcript.get_chrom()
        pos = str(var_dict["pos"] + 1)  # add 1 to make 1 based
        ref = str(var_dict["ref"])
        alt = str(var_dict["alt"])
        ID = "_".join(["str", pos, str(self.length)])
        if self.zygosity == "HOMOZYGOUS":
            zygosity = "1/1"
        if self.zygosity == "HEMIZYGOUS":
            zygosity = "1"
        if self.zygosity == "HETEROZYGOUS":
            zygosity = "0/1"
        return "\t".join([chrom, pos, ID, ref, alt, ".", ".", ".", "GT", zygosity])
