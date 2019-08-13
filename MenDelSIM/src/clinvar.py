from MenDelSIM.src.variant import Variant
from MenDelSIM.src.transcript import Transcript
import random

class ClinVar(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template: dict, transcript:Transcript):
        Variant.__init__(self, var_template, transcript)
        self.id = self.impact['clinvar_id']
        self.start = self.impact["START"] - 1
        self.ref = self.impact["REF"]
        self.alt = self.impact["ALT"]
        if self.type is not "ClinVar":
            raise ValueError("type must be ClinVar")
        if (len(self.ref) > 1):
            self.check_location(self, self.start, self.start + len(self.ref))
        elif (len(self.ref) == 1):
            self.check_location(self.start)

    def get_vcf_row(self) -> str:
        chrom = self.transcript.get_chr()
        pos = str(self.start + 1)  # add 1 to make 1 based
        ref = str(self.ref)
        alt = str(self.alt)
        ID = str(self.id)
        if self.zygosity == "HOMOZYGOUS":
            zygosity = "1/1"
        if self.zygosity == "HEMIZYGOUS":
            zygosity = "1"
        if self.zygosity == "HETEROZYGOUS":
            zygosity = "0/1"
        return "\t".join([chrom, pos, ID, ref, alt, ".", ".", ".", "GT", zygosity])
