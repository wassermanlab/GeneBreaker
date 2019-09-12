from MenDelSIM.src.variant import Variant
from MenDelSIM.src.transcript import Transcript
from lxml import etree


class CopyNumberVariant(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template: dict, transcript: Transcript):
        Variant.__init__(self, var_template, transcript)
        self.chrom = self.transcript.get_chr()
        self.start = self.impact["START"]
        self.end = self.impact["END"]
        self.copy_change = self.impact["COPY_CHANGE"]
        self.check_cnv()

    def check_cnv(self):
        """checks all cnv features and fixes to 0 based"""
        if self.type != "CNV":
            raise ValueError("Must be CNV type")
        if self.copy_change == 0 or self.copy_change < -1:
            raise ValueError(
                "Copy change cannot be less than -1 or equal to  0.")
        if type(self.start) != int or type(self.end) != int:
            raise ValueError("start and end must be int.")
        self.start = self.start-1
        if self.copy_change < 0:
            self.check_location(self.start, self.end)
        else:
            end = self.end + (self.end - self.start)*self.copy_change
            self.check_location(self.start, end)
        if (self.end - self.start) < 50:
            raise ValueError("cnvs must be at least 50 bp long")

    def get_anchor_position(self) -> str:
        """Retrieve a DNA sequence from UCSC.
        Note: UCSC assumes 1 based indexing so we add a 1"""
        # Initialize
        seq = self.get_seq(self.chrom, (self.start-1),
                           self.start, self.transcript.get_genome())
        return seq.upper()

    def get_region_seq(self) -> str:
        # Get sequence
        seq = self.get_seq(self.chrom, self.start, self.end,
                           self.transcript.get_genome())
        return seq.upper()

    def get_vcf_row(self, format: str = "simple") -> dict:
        # get regions
        chrom = self.chrom
        pos = str(self.start + 1)  # add 1 to make 1 based
        end = str(self.end)
        distance = int(end) - int(pos)
        if format == "simple":
            if self.copy_change > 0:
                ref = self.get_region_seq()
                alt = self.get_region_seq()*self.copy_change
            else:  # deletion
                pos = str(self.start)
                ref = self.get_anchor_position() + self.get_region_seq()
                alt = self.get_anchor_position()
            info = "."
        # TODO: fix 4.2 format
        # elif format == "4.2":
        #     ref = "N"
        #     if self.length > 0:  # duplication
        #         alt = "<DUP>"
        #         end = str(int(pos) + distance*self.copy_change)
        #         info = "SVTYPE=DUP;END=" + end + \
        #             ";SVLEN=-"+str(distance*self.copy_change)
        #     else:  # deletion
        #         alt = "<DEL>"
        #         info = "SVTYPE=DEL;END=" + end + ";SVLEN=-"+str(distance)
        ID = "_".join(["cnv", pos, str(self.copy_change)])
        if self.zygosity == "HOMOZYGOUS":
            zygosity = "1/1"
        if self.zygosity == "HEMIZYGOUS":
            zygosity = "1"
        if self.zygosity == "HETEROZYGOUS":
            zygosity = "0/1"

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
