from simulator.src.variant import Variant
from simulator.src.transcript import Transcript
import random
from GUD.ORM import ClinVar as CLN
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
from . import establish_GUD_session

class ClinVar(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template: dict):
        Variant.__init__(self, var_template)
        self.impact = self.impact["CLINVAR_UID"]
        print('check that is clinvar')

    def get_clinvar_information(self):
        """return the motif of the specific str"""
        session = establish_GUD_session()
        clinvar = CLN()
        clinvar = clinvar.select_by_name(session, self.impact, True)
        if clinvar is None:
            raise Exception("no variant with that clinvar ID")
        return {"pos": clinvar.start,
                "ref": clinvar.qualifiers["ref"],
                "alt": clinvar.qualifiers["alt"],
                "id": self.impact}

    def get_vcf_row(self, transcript: Transcript) -> str:
        # get regions
        # TODO: put check here that its in correct region??
        var_dict = self.get_clinvar_information()

        chrom = transcript.chrom
        pos = str(var_dict["pos"] + 1)  # add 1 to make 1 based
        ref = str(var_dict["ref"])
        alt = str(var_dict["alt"])
        ID = str(var_dict["id"])
        if self.zygosity == "HOMOZYGOUS":
            zygosity = "1/1"
        if self.zygosity == "HEMIZYGOUS":
            zygosity = "1"
        if self.zygosity == "HETEROZYGOUS":
            zygosity = "0/1"
        return "\t".join([chrom, pos, ID, ref, alt, ".", ".", ".", "GT", zygosity])
