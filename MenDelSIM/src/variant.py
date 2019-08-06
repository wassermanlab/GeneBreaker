import re
from MenDelSIM.src.transcript import Transcript
from lxml import etree
from MenDelSIM.src.api_helper import *



class Variant:
    # template is type dict
    def __init__(self, var_template: dict, transcript: Transcript):
        if (var_template is None):
            return None
        self.transcript = Transcript
        self.type = var_template["TYPE"]
        self.region = var_template["REGION"]
        self.impact = var_template["IMPACT"]
        self.zygosity = var_template["ZYGOSITY"]
        try: 
            self.check_variant()
        except Exception as e:
            raise(e)

    def parse_region(self) -> dict:
        """sets self.region to {chrom: ___, start: ___, stop: ___}"""
        chrom = var_region.split(":")[0]
        start = int(var_region.split(":")[1].split(
            "-")[0]) - 1  # making 0 based
        end = int(var_region.split(":")[1].split("-")[1])
        if start >= end:
            raise ValueError("start is greater than end in custom position")
        if start < 0:
            raise ValueError("invalid start position")
        if end > get_chromosome(self.transcript.chrom, self.transcript.genome)["size"]:
            raise ValueError(
                "invalid end position, greater than chromosome length")
        if chrom != transcript.chrom:
            raise ValueError(
                "custom location is not on the same chromosome as the transcript")
        return (chrom, (start, end))

    # TODO: fix occurrences of this
    def get_seq(self, chrom: int, start: int, end: int, genome: str) -> str:
        """Retrieve a DNA sequence from UCSC.
        Note: UCSC assumes 1 based indexing so we add a 0"""
        # Initialize
        sequence = ""
        url = "http://genome.ucsc.edu/cgi-bin/das/%s/dna?segment=%s:%s,%s" % (
            genome, chrom, start+1, end)
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
        if self.region not in ['CODING', 'UTR', 'INTRONIC', 'PROMOTER', 'ENHANCER', 'GENIC']:
            if re.match("^chr([XYM]|[1-9]|1[0-9]|2[0-2]):\d+-\d+$", self.region) is None:
                raise ValueError('region not one of the recognized regions: CODING, UTR, INTRONIC, PROMOTER, ENHANCER, GENIC, or the custom format')
            else:
                self.get_region()

    def check_type(self) -> bool: 
        if self.type not in ['SNV', 'INDEL', 'CNV', 'MEI', 'STR', 'ClinVar']:
            raise ValueError('TYPE not one of the recognized types: SNV, INDEL, CNV, MEI, STR, ClinVar')
        return True

    def check_zygosity(self): 
        if self.zygosity not in ['HOMOZYGOUS', 'HETEROZYGOUS', 'HEMIZYGOUS']:
            raise ValueError('ZYGOSITY must be one of the following: HOMOZYGOUS, HETEROZYGOUS, HEMIZYGOUS')
        return True 

    def check_variant(self): 
        try: 
            self.check_region()
            self.check_type()
            self.check_zygosity()
            return True 
        except Exception as e:
            raise(e)