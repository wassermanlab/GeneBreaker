import re
from MenDelSIM.src.transcript import Transcript
from lxml import etree


class Variant:
    # template is type dict
    def __init__(self, var_template: dict):
        self.type = var_template["TYPE"]
        self.region = var_template["REGION"]
        self.impact = var_template["IMPACT"]
        self.zygosity = var_template["ZYGOSITY"]
        
        if self.type not in ['SNV', 'INDEL', 'CNV', 'SV', 'MEI', 'STR', 'ClinVar']:
            raise ValueError('type not one of the recognized types: SNV, INDEL, CNV, SV, MEI, STR, ClinVar')
        
        if self.region not in ['CODING', 'UTR', 'INTRONIC', 'PROMOTER', 'ENHANCER', 'GENIC']:
            if re.match("^chr([XY]|[1-9]|1[0-9]|2[0-2]):\d+-\d+$", self.region) is None:
                raise ValueError('region not one of the recognized regions: CODING, UTR, INTRONIC, PROMOTER, ENHANCER, GENIC, or the custom format')
        
        if self.zygosity not in ['HOMOZYGOUS', 'HETEROZYGOUS', 'HEMIZYGOUS']:
            raise ValueError(
                'must have a zygosity within this list HOMOZYGOUS, HETEROZYGOUS, HEMIZYGOUS')

    def parse_region(transcript: Transcript, var_region: str) -> tuple:
        """returns a tuple (chrom, (start, stop))"""
        if re.match("^chr([XY]|[1-9]|1[0-9]|2[0-2]):\d+-\d+$", var_region) is None:
            raise Exception("custom region is not in correct format")
        chrom = var_region.split(":")[0]
        start = int(var_region.split(":")[1].split(
            "-")[0]) - 1  # making 0 based
        end = int(var_region.split(":")[1].split("-")[1])
        if start >= end:
            raise ValueError("start is greater than end in custom position")
        if start < 1:
            raise ValueError("invalid start position")
        if end > Chrom().chrom_size(session, chrom):
            raise ValueError(
                "invalid end position, greater than chromosome length")
        if chrom != transcript.chrom:
            raise ValueError(
                "custom location is not on the same chromosome as the transcript")
        return (chrom, (start, end))

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
