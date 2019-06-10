from MenDelSIM.src.variant import Variant
from MenDelSIM.src.transcript import Transcript 
import random
from Bio import SeqIO
import os
import zipfile

class MEI(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template: dict):
        Variant.__init__(self, var_template)
        self.element = self.impact["ELEMENT"]
        self.location = self.impact["LOCATION"]
        self.element_dict = {"ALU_MELT": "ALU", "LINE1_MELT": "LINE1", "SVA_MELT": "SVA"}
        if self.location != "ANY" and type(self.location) is not int: 
            raise ValueError("location must be ANY or int")
        if self.element not in ["ALU_MELT", "LINE1_MELT", "SVA_MELT"]:
            raise ValueError("""Only elements in the set of, "ALU_MELT", "LINE1_MELT", "SVA_MELT"
            are currently supported.""")
        if self.type != "MEI":
            raise ValueError("Must be MEI type")


    def get_insertion_str(self) -> str:
        """reads the fasta and gets string of fasta file"""
        ## get file name
        THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
        THIS_FOLDER = os.path.split(THIS_FOLDER)[0]
        zip_file = THIS_FOLDER + '/static/'+ self.element + ".zip"
        fasta_file_full = THIS_FOLDER + '/static/' + self.element_dict[self.element] + ".fa"
        fasta_file = self.element_dict[self.element]+".fa"
        ## unzip file 
        zip_file = zipfile.ZipFile(zip_file)
        zip_file.extract(fasta_file, os.path.join(THIS_FOLDER, 'static'))
        zip_file.close()
        
        insertion = ""
        for seq_record in SeqIO.parse(fasta_file_full, "fasta"):
            insertion = str(seq_record.seq)

        os.remove(fasta_file_full)
        return insertion.upper() 

    def get_insertion(self, transcript: Transcript) -> dict:
        """returns (ref, alt) tuple of insersion"""
        # get requested region
        if self.region in ["CODING", "INTRONIC", "UTR", "GENIC"]:
            regions = transcript.get_requested_region(self.region)
            if len(regions) == 0:
                raise ValueError("region requested does not exists")
        else: 
            regions = [self.parse_region(transcript, self.region)[1]]
        # get ranges
        region_range = []
        for region in regions:
            region_range = region_range + list(range(region[0], region[1]))
        if self.location == "ANY": # pick any position within the ranges
            pos = random.choice(region_range)
        else:
            if self.location not in region_range:
                raise ValueError("position must be within range")
            pos = self.location
        return {"pos": pos,
                "ref": self.get_seq(transcript.chrom, pos, pos+1),
                "alt": self.get_seq(transcript.chrom, pos, pos+1) + self.get_insertion_str()}

    def get_vcf_row(self, transcript: Transcript) -> str:
        chrom = str(transcript.get_chr())
        var_dict = self.get_insertion(transcript)
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
