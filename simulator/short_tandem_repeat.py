from simulator.variant import Variant
from simulator.transcript import Transcript 
import random
from GUD2.ORM import ShortTandemRepeat as GUDSTR ##
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session

class ShortTandemRepeat(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template):
        try:
            Variant.__init__(self, var_template)
            self.chrom = self.impact["CHROM"] 
            self.start = self.impact["START"]
            self.end = self.impact["END"]
            self.length = self.impact["STR"]
            if self.length == 0:
                raise Exception("STR length must be not equal to 0")
            if self.type != "STR":
                raise Exception("Must be STR type")
        except:
            print('check that type is STR and that impact is correctly formatted')
    
    def get_str_motif(self):
        """return the motif of the specific str"""
        db_name = "mysql://{}:@{}:{}/{}".format("ontarget_r", #todo change this to hg19
                                        "ontarget.cmmt.ubc.ca", "5506", "tamar_test")
        engine = create_engine(db_name, echo=False)
        session = Session(engine) 
        STR = GUDSTR()
        STR = STR.select_by_exact_location(session, self.chrom, self.start, self.end)
        return STR[0].motif

    def get_retraction(self):
        """returns (pos ,ref, alt) tuple of retraction"""
        #get str
        STR = self.get_str_motif()
        size = len(STR)
        total_repeat_length = self.end - self.start 
        #check that requested size is not over
        if total_repeat_length < -1 *size*self.length:
            raise Exception("retraction length is larger than the total str")
        positive_len = (self.length*-1)
        ref = STR*positive_len
        alt = ""
        return {"pos": self.start,
                "ref": ref,
                "alt": alt}

    def get_expantion(self):
        """returns (pos ,ref, alt) tuple of retraction"""
        #get str
        STR = self.get_str_motif()
        size = len(STR)
        if size*self.length> 20000:
            raise Exception("expansion length is too large")
        ref = STR[0]
        alt = STR*self.length + STR[0]
        return {"pos": self.start,
                "ref": ref,
                "alt": alt}

    def get_vcf_row(self, transcript):
        # get regions 
        regions = transcript.get_requested_region(self.region)
        check = False
        for region in regions:
            if region[0] <= self.start < region[1]:
                check = True
        if check is False: 
            raise Exception("STR is not in requested region")
        if self.length > 0:  # insersion
            var_dict = self.get_expantion()
        if self.length < 0: # deletion
            var_dict = self.get_retraction()
        pos = str(var_dict["pos"] + 1) # add 1 to make 1 based
        ref = str(var_dict["ref"])
        alt = str(var_dict["alt"])
        ID = "_".join(["str", pos, str(self.impact)])
        return "\t".join([chrom, pos, ID, ref, alt])
