from GeneBreaker.src.variant import Variant
from GeneBreaker.src.transcript import Transcript
import random
from GeneBreaker.src.api_helper import *
from lxml import etree


class ShortTandemRepeat(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template, transcript: Transcript):
        try:
            Variant.__init__(self, var_template, transcript)
            self.chrom = self.transcript.get_chr()
            STR = get_str(self.transcript.get_genome(), int(self.impact["STR_ID"]))
            self.start = int(STR["start"])
            self.end = int(STR["end"])
            self.length = self.impact["LENGTH"]
            self.motif = STR["qualifiers"]["motif"]
            self.check_str()
            # get deletion/insersion 
            if self.length > 0:  # insersion
                var_dict = self.get_expantion()
            if self.length < 0:  # deletion
                var_dict = self.get_retraction()
            self.pos = var_dict["pos"]
            self.ref = str(var_dict["ref"])
            self.alt = str(var_dict["alt"])
            self.id = "_".join(["str", str(self.pos), str(self.length)])
        except Exception as e:
            raise(e)
    
    def check_str(self):
        if self.type != "STR":
            raise ValueError("Must be STR type")
        if type(self.start) != int or type(self.end) != int:
            raise ValueError("start and end must be int.")
        if (self.length>0):
            self.check_location(self.start)
        else: 
            end = self.start + self.length*-1*len(self.motif)
            self.check_location(self.start-1,end)
        
    def get_anchor_position(self) -> str:
        """Retrieve a DNA sequence from UCSC.
        Note: UCSC assumes 1 based indexing so we add a 1"""
        seq = self.get_seq(self.chrom, (self.start-1), self.start, self.transcript.get_genome())
        return seq.upper()

    def get_retraction(self) -> dict:
        """returns (pos ,ref, alt) tuple of retraction"""
        # get str
        size = len(self.motif)
        total_repeat_length = self.end - self.start
        # check that requested size is not over
        positive_len = (self.length*-1)
        if total_repeat_length < positive_len*size:
            raise Exception("retraction length is larger than the total str")
        alt = self.get_anchor_position()
        ref = alt + self.motif*positive_len
        return {"pos": self.start-1,
                "ref": ref,
                "alt": alt}

    def get_expantion(self) -> dict:
        """returns (pos ,ref, alt) tuple of retraction"""
        size = len(self.motif)
        if size*self.length > 20000:
            raise Exception("expansion length is too large")
        ref = self.motif
        alt = self.motif*self.length + self.motif
        return {"pos": self.start,
                "ref": ref,
                "alt": alt}
