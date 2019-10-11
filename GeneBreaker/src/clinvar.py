from GeneBreaker.src.variant import Variant
from GeneBreaker.src.transcript import Transcript
import random
from GeneBreaker.src.api_helper import *


## TODO: remove ref and alt from here
class ClinVar(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template: dict, transcript:Transcript):
        Variant.__init__(self, var_template, transcript)
        self.id = self.impact["CLINVAR_ID"]
        clinvar = get_clinvar(self.transcript.get_genome(), self.id)
        if clinvar != []:
            self.pos = clinvar['start']
            self.ref = clinvar['qualifiers']['ref']
            self.alt = clinvar['qualifiers']['alt']
        else: ## checks
            raise ValueError("type must be CLINVARID_invalid") #no result
        if self.type != "CLINVAR":
            raise ValueError("type must be CLINVAR")
        ## bad location
        if (len(self.ref) > 1):
            self.check_location(self.pos, self.pos + len(self.ref))
        elif (len(self.ref) == 1):
            self.check_location(self.pos)
