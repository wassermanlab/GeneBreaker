import json
import argparse
import json
from GUD.ORM import Gene
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
from lxml import etree

class Transcript:

    def __init__(self, name, chrom, txStart, txEnd, strand):
        """ create transcript object 
        """
        # Establish a SQLalchemy session w/ GUD
        db_name = "mysql://{}:@{}:{}/{}".format("ontarget_r",
            "ontarget.cmmt.ubc.ca", "5506", "hg19")

        try:
            self.name = name
            self.chrom = chrom
            self.txStart = txStart
            self.txEnd = txEnd
            self.strand = strand

            engine = create_engine(db_name, echo=False)
            session = Session(engine)
            gene = Gene()
            transcripts = Gene.select_by_name(session, self.name)
            for t in transcripts:
                if str(t) == "{}\t{}\t{}\t{}\t0\t{}".format(self.chrom, self.txStart, self.txEnd, self.name, self.strand):
                    gene = t
            self.cdsStart = gene.cdsStart
            self.cdsEnd = gene.cdsEnd
            self.exonStarts = gene.exonStarts
            self.exonEnds = gene.exonEnds
        except:
            print("Unexpected error, check validity of config and gene inputs")

    def get_seq(self): 
        """Retrieve a DNA sequence from UCSC.
        Note: UCSC assumes 1 based indexing so we add a 0""" 
        # Initialize
        sequence = ""
        url = "http://genome.ucsc.edu/cgi-bin/das/%s/dna?segment=%s:%s,%s" % (
            "hg19", self.chrom, self.txStart+1, self.txEnd)
        # Get XML
        xml = etree.parse(url, parser=etree.XMLParser())
        # Get sequence
        sequence = xml.xpath("SEQUENCE/DNA/text()")[0].replace("\n", "")
        return sequence.upper()

    def get_chr(self): 
        """ returns chr from gene""" 
        return self.chrom       

    def get_exons(self):
        """ returns the codons from a gene [(start, stop)*] sorted from lowest to highest start position"""
        starts = self.exonStarts.split(",")[:-1]
        starts.sort(key=int)
        ends = self.exonEnds.split(",")[:-1]
        ends.sort(key=int)
        return zip(starts ,ends)
        

    def get_coding(self):
        """ returns all coding regions from gene in the form of [[start, stop],*]"""
        exons = self.get_exons()
        coding_exons = []
        for exon in exons: 
            if exon[0] > self.cdsStart and exon[1]>

    def get_introns(self):
        """ returns all introns from gene [[start, stop],*]"""

    def get_utr(self):
        """ returns the utrs """

    def get_requested_region(self, region):
        if region == "CODING":
            return self.get_coding()
        elif region == "UTR":
            return self.get_utr()
        elif region == "INTRONIC":
            return self.get_introns()    
        else:
            raise Exception("region not valid")

    def get_codon_from_pos(self, pos):
        """ from position get codon matching to position
        return codon and position of nucleotide in codon (codon, pos) """
        # Todo: code entire logic and test
        coding_regions = self.get_coding()
        coding_regions = sorted(coding_regions, key=lambda x: x[0])
        gene = ""
        shift_pos = 0
        last_exon = None
        coding_pos = None       # will be used to get the position requested without introns
        for exon in coding_regions:
            gene = gene + self.seq[exon[0]:exon[1]]
            if last_exon is None:
                last_exon = exon[1]
            else:
                shift_pos = shift_pos + exon[0] - last_exon
                last_exon = exon[1]
            if exon[0] <= pos < exon[1]:
                coding_pos = pos - shift_pos
        if coding_pos is None:
            raise Exception("position specified is not in coding region.")
        position_in_codon = coding_pos % 3
        if coding_pos < 3:
            return(self.seq[0:3], position_in_codon)
        else:
            start = coding_pos - position_in_codon
            stop = start + 3
            return(gene[start:stop], position_in_codon)

    def __str__(self):
        return "{}\t{}\t{}\t{}\t0\t{}".format(self.chrom, self.txStart, self.txEnd, self.name, self.strand)