import re
from GeneBreaker.src.transcript import Transcript
from lxml import etree
from GeneBreaker.src.api_helper import *


class Variant:
    # template is type dict
    def __init__(self, var_template: dict, transcript: Transcript):
        try:
            if (var_template == None):
                return None
            self.transcript = transcript
            self.type = var_template["TYPE"]
            self.region = var_template["REGION"]
            self.impact = var_template["IMPACT"]

            self.chrom = self.transcript.get_chr()
            self.pos = None
            self.id = None
            self.ref = None
            self.alt = None
            self.zygosity = var_template["ZYGOSITY"]
            self.check_variant()
        except Exception as e:
            raise(e)

    def check_location(self, start, end=None):
        # assume 0 based at this point
        """checks location for validity"""
        region_map = self.transcript.get_regionmap()
        passing = False
        if end == None:
            for i in region_map[self.region]:
                if start <= i[1] and start >= i[0]:
                    passing = True
        else:
            if self.region == "GENIC" or self.region == "CODING":
                for i in region_map[self.region]:
                    if start < i[1] and end > i[0]:
                        passing = True
            elif self.region == "UTR":
                for i in region_map["UTR"]:
                    if start < i[1] and end > i[0]:
                        passing = True
                for i in region_map["CODING"]:
                    if start < i[1] and end > i[0]:
                        passing = False
            elif self.region == "INTRONIC":
                for i in region_map["INTRONIC"]:
                    if start >= i[0] and end <= i[1]:
                        passing = True
            else:
                i = self.parse_region()
                if start >= i[1][0] and end <= i[1][1]:
                    passing = True

        if not passing:
            raise ValueError("invalid position.")

    def parse_region(self) -> dict:
        """sets self.region to {chrom: ___, start: ___, stop: ___}"""
        chrom = self.region.split(":")[0]
        start = int(self.region.split(":")[1].split(
            "-")[0]) - 1  # making 0 based
        end = int(self.region.split(":")[1].split("-")[1])
        if start >= end:
            raise ValueError("start is greater than end in custom position")
        if start < 0:
            raise ValueError("invalid start position")
        if end > get_chromosome(self.transcript.chrom, self.transcript.genome)["size"]:
            raise ValueError(
                "invalid end position, greater than chromosome length")
        if chrom != self.transcript.chrom:
            raise ValueError(
                "custom location is not on the same chromosome as the transcript")
        return (chrom, (start, end))

    def get_seq(self, chrom: int, start: int, end: int, genome: str) -> str:
        """Retrieve a DNA sequence from UCSC.
        Note: UCSC assumes 1 based indexing so we add a 1"""
        # Initialize
        sequence = ""
        if (genome == "grch38"):
            g = "hg38"
        elif (genome == "grch37"):
            g = "hg19"
        url = "http://genome.ucsc.edu/cgi-bin/das/%s/dna?segment=%s:%s,%s" % (
            g, chrom, start+1, end)
        # Get XML
        xml = etree.parse(url, parser=etree.XMLParser())
        # Get sequence
        sequence = xml.xpath("SEQUENCE/DNA/text()")[0].replace("\n", "")
        sequence = str(sequence)
        return sequence.upper()

    def get_type(self) -> str:
        return self.type

    def get_region(self) -> str:
        return self.region

    def check_region(self) -> bool:
        """checks region and gets custom region"""
        if self.region not in ['CODING', 'UTR', 'INTRONIC', 'PROMOTER', 'ENHANCER', 'GENIC']:
            if re.match("^([XYM]|[1-9]|1[0-9]|2[0-2]):\d+-\d+$", self.region) == None:
                raise ValueError(
                    'region not one of the recognized regions: CODING, UTR, INTRONIC, PROMOTER, ENHANCER, GENIC, or the custom format')
            else:
                self.get_region()

    def check_type(self) -> bool:
        """checks type"""
        if self.type not in ['SNV', 'INDEL', 'CNV', 'MEI', 'STR', 'CLINVAR']:
            raise ValueError(
                'TYPE not one of the recognized types: SNV, INDEL, CNV, MEI, STR, CLINVAR')
        return True

    def check_zygosity(self):
        """checks that zygosity is valid"""
        if self.zygosity not in ['HOMOZYGOUS', 'HETEROZYGOUS', 'HEMIZYGOUS']:
            raise ValueError(
                'ZYGOSITY must be one of the following: HOMOZYGOUS, HETEROZYGOUS, HEMIZYGOUS')
        return True

    def check_variant(self):
        """checks variant for validity"""
        try:
            self.check_region()
            self.check_type()
            self.check_zygosity()
            return True
        except Exception as e:
            raise(e)

    def get_region_range(self):
        """get ther region requested for that variant using the attached transcript"""
        regions = self.transcript.get_requested_region(self.region)
        region_range = []
        for region in regions:
            region_range = region_range + list(range(region[0], region[1]))
        if len(region_range) == 0:
            raise ValueError("""invalid region selected""")
        return region_range

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
        return {
            "chrom": chrom,
            "pos":  pos,
            "id": ID,
            "ref": ref,
            "alt": alt,
            "qual": ".",
            "filter": ".",
            "info": ".",
            "format": "GT",
            "proband": zygosity}
