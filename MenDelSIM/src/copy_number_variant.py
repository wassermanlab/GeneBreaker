from MenDelSIM.src.variant import Variant
from MenDelSIM.src.transcript import Transcript
from lxml import etree


class CopyNumberVariant(Variant):
    # assume var_template is of type dict already
    def __init__(self, var_template: dict, transcript: Transcript):
        Variant.__init__(self, var_template, transcript)
        self.chrom = self.impact["CHROM"]
        self.start = self.impact["START"]
        self.end = self.impact["END"]
        self.length = self.impact["COPY_CHANGE"]
        self.check_cnv()



    # TODO: check location
    def check_location(self):
        """checks location and converts start to 0 based"""
        self.start = self.start - 1 
        self.transcript.get_regionmap()


    def check_cnv(self):
        """checks all cnv features and fixes to 0 based"""
        if self.type != "CNV":
            raise ValueError("Must be CNV type")
        self.check_location()
        if self.length == 0 or self.length < -1 or self.length == 1:
            raise Exception("Copy change cannot be less than -1 or equal to 1 or 0.")

    def get_anchor_position(self) -> str:
        """Retrieve a DNA sequence from UCSC.
        Note: UCSC assumes 1 based indexing so we add a 1"""
        # Initialize
        seq = self.get_seq(self.chrom, (self.start-1), self.start, self.transcript.get_genome())
        return seq.upper()

    def get_region_seq(self) -> str:
        """Retrieve a DNA sequence from UCSC.
        Note: UCSC assumes 1 based indexing so we add a 1"""
        # Get sequence
        seq = self.get_seq(self.chrom, self.start, self.end, self.transcript.get_genome())
        return seq.upper()

    def get_vcf_row(self, transcript: Transcript, format: str = "simple") -> str:
        # get regions
        chrom = self.chrom
        pos = str(self.start + 1)  # add 1 to make 1 based
        end = str(self.end)
        distance = int(end) - int(pos)
        if format == "simple":
            if self.length > 0:
                ref = self.get_region_seq()
                alt = self.get_region_seq()*self.length
            else:  # deletion
                ref = self.get_anchor_position() + self.get_region_seq()
                alt = self.get_anchor_position()
            info = "."
        elif format == "4.2":
            ref = "N"
            if self.length > 0:  # duplication
                alt = "<DUP>"
                end = str(int(pos) + distance*self.length)
                info = "SVTYPE=DUP;END=" + end + \
                    ";SVLEN=-"+str(distance*self.length)
            else:  # deletion   
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
