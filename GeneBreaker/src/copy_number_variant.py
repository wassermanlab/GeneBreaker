from GeneBreaker.src.variant import Variant
from GeneBreaker.src.transcript import Transcript
from lxml import etree


class CopyNumberVariant(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template: dict, transcript: Transcript):
        try:
            Variant.__init__(self, var_template, transcript)
            self.chrom = self.transcript.get_chr()
            self.pos = self.impact["START"]
            self.end = self.impact["END"]
            self.copy_change = self.impact["COPY_CHANGE"]
            self.check_cnv()
            self.ref = self.get_anchor_position()
            if self.copy_change < 0 : # deletion
                self.alt = "<DEL>"
            else: 
                self.alt = "<DUP:TANDEM>"
            self.id = "CNV_" + str(self.pos) 
        except Exception as e:
            raise(e)

    def check_cnv(self):
        """checks all cnv features and fixes to 0 based"""
        if self.type != "CNV":
            raise ValueError("Must be CNV type")
        if self.copy_change not in  [-1, 1]:
            raise ValueError(
                "Copy change must be -1 or 1")
        if type(self.pos) != int or type(self.end) != int:
            raise ValueError("start and end must be int.")
        self.pos = self.pos-1
        if self.copy_change < 0:
            self.check_location(self.pos, self.end)
        else:
            end = self.end + (self.end - self.pos)*self.copy_change
            self.check_location(self.pos, end)
        if (self.end - self.pos) < 50:
            raise ValueError("cnvs must be at least 50 bp long")

    def get_anchor_position(self) -> str:
        """Retrieve a DNA sequence from UCSC.
        Note: UCSC assumes 1 based indexing so we add a 1"""
        # Initialize
        seq = self.get_seq(self.chrom, self.pos,
                           self.pos+1, self.transcript.get_genome())
        return seq.upper()

    def get_vcf_row(self) -> dict:
        chrom = self.chrom
        pos = str(self.pos + 1)  # add 1 to make 1 based
        ref = str(self.ref)
        alt = str(self.alt)
        ID = str(self.id)
        if self.zygosity == "HOMOZYGOUS":
            zygosity = "1/1"
        if self.zygosity == "HEMIZYGOUS":
            zygosity = "1/1"
        if self.zygosity == "HETEROZYGOUS":
            zygosity = "0/1"
        svend = self.end
        svlen = self.end - self.pos + 1 
        if self.copy_change < 0: # deletion 
            info = "SVTYPE=DEL;END=" + str(svend) + ";SVLEN=-" + str(svlen) + ";"
        else: 
            info = "SVTYPE=DUP;END=" + str(svend) + ";SVLEN=" + str(svlen) + ";"
        return {
            "chrom": chrom,
            "pos":  pos,
            "id": ID,
            "ref": ref,
            "alt": alt,
            "qual": ".",
            "filter": ".",
            "info": info,
            "format": "GT",
            "proband": zygosity}
