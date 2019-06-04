from simulator.variant import Variant
from simulator.transcript import Transcript 
import random
from GUD2.ORM import CNV
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
from lxml import etree

class CopyNumberVariant(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template):
        Variant.__init__(self, var_template)
        self.chrom = self.impact["CHROM"] 
        self.start = self.impact["START"]
        self.end = self.impact["END"]
        self.length = self.impact["CNV"]
        if self.length == 0:
            raise Exception("CNV length must be not equal to 0")
        if self.type != "CNV":
            raise Exception("Must be CNV type")
    
    def get_anchor_position(self):
        """Retrieve a DNA sequence from UCSC.
        Note: UCSC assumes 1 based indexing so we add a 1""" 
        # Initialize
        sequence = ""
        url = "http://genome.ucsc.edu/cgi-bin/das/%s/dna?segment=%s:%s,%s" % (
            "hg19", self.chrom, self.start, self.start) ## by using the start without adjusting we are getting the position prior to the start
        # Get XML
        xml = etree.parse(url, parser=etree.XMLParser())
        # Get sequence
        sequence = xml.xpath("SEQUENCE/DNA/text()")[0].replace("\n", "")
        return sequence.upper()

    def get_region_seq(self):
        """Retrieve a DNA sequence from UCSC.
        Note: UCSC assumes 1 based indexing so we add a 1""" 
        # Initialize
        sequence = ""
        url = "http://genome.ucsc.edu/cgi-bin/das/%s/dna?segment=%s:%s,%s" % (
            "hg19", self.chrom, self.start+1, self.end) ## by using the start without adjusting we are getting the position prior to the start
        # Get XML
        xml = etree.parse(url, parser=etree.XMLParser())
        # Get sequence
        sequence = xml.xpath("SEQUENCE/DNA/text()")[0].replace("\n", "")
        return sequence.upper()

    def get_vcf_row(self, transcript, format = "simple"):
        # get regions 
        chrom = self.chrom
        pos = str(self.start + 1) # add 1 to make 1 based
        end = str(self.end)
        distance = int(end) - int(pos)
        if format == "simple":
            if self.length > 0:
                ref = self.get_region_seq()
                alt = self.get_region_seq()*self.length
            else: ## deletion
                ref = self.get_anchor_position() + self.get_region_seq()
                alt = self.get_anchor_position()
            info = "."
        elif format == "4.2":
            ref = "N"
            if self.length > 0: ## duplication
                alt = "<DUP>"
                end = str(int(pos) + distance*self.length)
                info = "SVTYPE=DUP;END=" + end + ";SVLEN=-"+str(distance*self.length)
            else: ## deletion   
                alt = "<DEL>"
                info = "SVTYPE=DEL;END=" + end + ";SVLEN=-"+str(distance)
        ID = "_".join(["cnv", pos, str(self.length)])
        if self.zygosity == "HOMOZYGOUS":
            zygosity = "1/1"  
        if self.zygosity == "HEMIZYGOUS":
            zygosity = "1"
        if self.zygosity == "HETEROZYGOUS":
            zygosity = "0/1"  
        
        return "\t".join([chrom, pos, ID, ref, alt, ".", ".", info, "GT", zygosity])
