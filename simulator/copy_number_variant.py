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
    
    def get_vcf_row(self, transcript):
        # get regions 
        chrom = self.chrom
        pos = str(self.start + 1) # add 1 to make 1 based
        end = str(self.end)
        distance = int(end) - int(pos)
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
