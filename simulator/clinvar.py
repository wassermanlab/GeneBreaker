from simulator.variant import Variant
from simulator.transcript import Transcript
import random
from GUD2.ORM import ClinVar as CLN
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session


class ClinVar(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template):
        try:
            Variant.__init__(self, var_template)
        except:
            print('check that is clinvar')

    def get_clinvar_information(self):
        """return the motif of the specific str"""
        db_name = "mysql://{}:@{}:{}/{}".format("ontarget_r",  # todo change this to hg19
                                                "ontarget.cmmt.ubc.ca", "5506", "tamar_test")
        engine = create_engine(db_name, echo=False)
        session = Session(engine)
        clinvar = CLN()
        clinvar = clinvar.select_by_name(session, self.impact)
        if clinvar is None:
            raise Exception("no variant with that clinvar ID")
        return {"pos": clinvar[1].start,
                "ref": clinvar[0].ref,
                "alt": clinvar[0].alt,
                "id": self.impact}

    def get_vcf_row(self, transcript):
        # get regions
        ## TODO: put check here that its in correct region??
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
