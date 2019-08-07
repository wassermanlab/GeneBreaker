from MenDelSIM.src.variant import Variant
from MenDelSIM.src.transcript import Transcript 
import random
from Bio import SeqIO
import os
from zipfile import ZipFile

class MEI(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template: dict, transcript: Transcript):
        Variant.__init__(self, var_template, transcript)
        self.element = self.impact["ELEMENT"]
        self.location = self.impact["LOCATION"]
        try: 
            self.check_mei()
        except Exception as e:
            raise(e)

    def check_element(self):
        """checks element validity"""
        if self.element not in ["ALU", "LINE", "SVA"]:
            raise ValueError("""Only elements in the set of: "ALU", "LINE", "SVA" are currently supported.""")

    def check_mei(self):
        """checks all mei features"""
        if self.type != "MEI":
            raise ValueError("Must be MEI type")
        self.check_element()
        self.check_location(self.location)

    def get_insertion_str(self) -> str:
        """reads the fasta and gets string of fasta file"""
        ## get file name
        THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
        THIS_FOLDER = os.path.split(THIS_FOLDER)[0]
        zip_file = THIS_FOLDER + '/static/'+ self.element + ".zip"
        fasta_file_full = THIS_FOLDER + '/static/' + self.element + ".fa"
        fasta_file = self.element +".fa"
        ## unzip file 
        with ZipFile(zip_file,"r") as zip_ref:
            zip_ref.extract(fasta_file, THIS_FOLDER + '/static/')
        
        insertion = ""
        for seq_record in SeqIO.parse(fasta_file_full, "fasta"):
            insertion = str(seq_record.seq)

        os.remove(fasta_file_full)
        return insertion.upper() 

    def get_insertion(self) -> dict:
        """returns (ref, alt) tuple of insersion in 0 based"""
        # get requested region range 
        region_range = self.get_region_range()
        if self.location == "ANY": # pick any position within the ranges
            pos = random.choice(region_range)
        else:
            pos = self.location
        return {"pos": pos,
                "ref": self.get_seq(self.transcript.chrom, pos, pos+1, self.transcript.genome),
                "alt": self.get_seq(self.transcript.chrom, pos, pos+1, self.transcript.genome) + self.get_insertion_str()}

    def get_vcf_row(self) -> str:
        chrom = str(self.transcript.get_chr())
        var_dict = self.get_insertion()
        pos = str(var_dict["pos"] + 1) # add 1 to make one based
        ref = str(var_dict["ref"])
        alt = str(var_dict["alt"])
        ID = "_".join(["mei", pos, str(self.element)])
        if self.zygosity == "HOMOZYGOUS":
            zygosity = "1/1"  
        if self.zygosity == "HEMIZYGOUS":
            zygosity = "1"
        if self.zygosity == "HETEROZYGOUS":
            zygosity = "0/1"  
        return "\t".join([chrom, pos, ID, ref, alt, ".", ".", ".", "GT", zygosity])
