from MenDelSIM.src.variant import Variant
from MenDelSIM.src.transcript import Transcript 
import random
from GUD.ORM import ShortTandemRepeat as GUDSTR ##
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
from lxml import etree
from . import establish_GUD_session

class ShortTandemRepeat(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template):
        Variant.__init__(self, var_template)
        self.chrom = self.impact["CHROM"] 
        self.start = self.impact["START"]
        self.end = self.impact["END"]
        self.length = self.impact["STR"]
        if self.length == 0:
            raise Exception("STR length must be not equal to 0")
        if self.type != "STR":
            raise Exception("Must be STR type")
    
    def get_str_motif(self) -> str:
        """return the motif of the specific str"""
        session = establish_GUD_session()
        STR = GUDSTR()
        STR = STR.select_by_exact_location(session, self.chrom, self.start, self.end, True)
        return STR.qualifiers["motif"]


    def get_anchor_position(self, chrom: str) -> str:
        """Retrieve a DNA sequence from UCSC.
        Note: UCSC assumes 1 based indexing so we add a 1""" 
        # Initialize
        sequence = ""
        url = "http://genome.ucsc.edu/cgi-bin/das/%s/dna?segment=%s:%s,%s" % (
            "hg19", chrom, self.start, self.start) ## by using the start without adjusting we are getting the position prior to the start
        # Get XML
        xml = etree.parse(url, parser=etree.XMLParser())
        # Get sequence
        sequence = xml.xpath("SEQUENCE/DNA/text()")[0].replace("\n", "")
        return sequence.upper()

    def get_retraction(self, chrom) -> dict:
        """returns (pos ,ref, alt) tuple of retraction"""
        #get str
        STR = self.get_str_motif()
        size = len(STR)
        total_repeat_length = self.end - self.start 
        #check that requested size is not over
        if total_repeat_length < -1 *size*self.length:
            raise Exception("retraction length is larger than the total str")
        positive_len = (self.length*-1)
        alt = self.get_anchor_position(chrom) 
        ref = alt + STR*positive_len
        return {"pos": self.start-1,
                "ref": ref,
                "alt": alt}

    def get_expantion(self) -> dict:
        """returns (pos ,ref, alt) tuple of retraction"""
        STR = self.get_str_motif()
        size = len(STR)
        if size*self.length> 20000:
            raise Exception("expansion length is too large")
        ref = STR[0]
        alt = STR*self.length + STR[0]
        return {"pos": self.start,
                "ref": ref,
                "alt": alt}

    def get_vcf_row(self, transcript: Transcript) -> str:
        # get regions 
        regions = transcript.get_requested_region(self.region) 
        check = False
        for region in regions:
            if region[0] <= self.start <= region[1]:
                check = True
        if check is False: 
            raise Exception("STR is not in requested region")
        if self.length > 0:  # insersion
            var_dict = self.get_expantion()
        if self.length < 0: # deletion
            var_dict = self.get_retraction(transcript.chrom)
        chrom = transcript.chrom
        pos = str(var_dict["pos"] + 1) # add 1 to make 1 based
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
