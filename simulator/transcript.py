import json
import argparse
import json
from GUD2.ORM import Gene
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
from lxml import etree
from Bio.Seq import Seq

class Transcript:

    def __init__(self, uid):
        """ create transcript object 
        """
        # Establish a SQLalchemy session w/ GUD
        db_name = "mysql://{}:@{}:{}/{}".format("ontarget_r",
            "ontarget.cmmt.ubc.ca", "5506", "tamar_test") #todo change back to hg19 once inserted

        try:
            engine = create_engine(db_name, echo=False)
            session = Session(engine)
            gene = Gene()
            transcript_region = Gene.select_by_uid_joined(session, uid)
            transcript = transcript_region[0]
            region = transcript_region[1]
            self.name = transcript.uid
            self.cdsStart   = long(transcript.cdsStart)
            self.cdsEnd     = long(transcript.cdsEnd)
            self.exonStarts = transcript.exonStarts
            self.exonEnds   = transcript.exonEnds
            self.chrom      = region.chrom
            self.txStart    = long(region.start)
            self.txEnd      = long(region.end)
            self.strand     = transcript.strand
        except:
            print("Unexpected error, check validity of config and gene inputs")

    def get_seq(self, stranded=False): 
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
        if stranded is False: 
            return sequence.upper()
        elif self.strand == "+":
            return sequence.upper()
        elif self.strand == "-":
            sequence = Seq(sequence).reverse_complement()
            sequence = str(sequence)
            return sequence.upper()

    def get_seq_from_pos(self, positions):
        """returns positive strand positions"""
        seq = self.get_seq()
        cut = ""
        if type(positions) is tuple:
            start = positions[0]-self.txStart
            end = positions[1]-self.txStart
            cut =  seq[start:end]
        elif type(positions) is list:
            for pos in positions:
                start = pos[0]-self.txStart
                end = pos[1]-self.txStart
                cut = cut + seq[start:end]
        else:
            raise Exception("wrong type to get_seq_from_pos(), should only except tuple or list or tuples")
        if self.strand == "+":
            return cut
        elif self.strand == "-":
            return str(Seq(cut).reverse_complement())

    def get_start(self):
        return self.txStart

    def get_chr(self): 
        """ returns chr from gene""" 
        return self.chrom       

    def positive_sorted(self, positions):
        positions.sort(key=lambda tup: tup[0], reverse=False)
        return positions 

    def negative_sorted(self, positions):
        positions.sort(key=lambda tup: tup[0], reverse=True)
        return positions 
   
    def get_exons(self):
        """ returns the codons from a gene [(start, stop)*] sorted from lowest to highest start position"""
        starts = self.exonStarts.split(",")[:-1]
        starts = [long(x) for x in starts]
        starts.sort()
        ends = self.exonEnds.split(",")[:-1]
        ends = [long(x) for x in ends]
        ends.sort()

        return zip(starts, ends) 
        
    def get_coding(self):
        """ returns all coding regions from gene in the form of [(start, stop),*]"""
        exons = self.get_exons()
        coding_exons = []
        if self.cdsStart == self.cdsEnd == self.txEnd:
            return coding_exons
        for exon in exons: 
            if exon[0] >= self.cdsStart and exon[1] <= self.cdsEnd: # exon is fully within coding region
                coding_exons.append(exon)
            elif exon[0] < self.cdsStart and self.cdsStart < exon[1] and exon[1] <= self.cdsEnd: #exon start is out of range but exon end is withing range
                coding_exons.append((self.cdsStart, exon[1]))
            elif exon[0] >= self.cdsStart and exon[0] < self.cdsEnd and  exon[1] > self.cdsEnd: #exon end is out of range but start is within range
                coding_exons.append((exon[0], self.cdsEnd))
            elif exon[0] < self.cdsStart and exon[1] > self.cdsEnd: #exon end and start are out of range 
                coding_exons.append((self.cdsStart, self.cdsEnd))
            else: #exon is completely out of range
                None 
        coding_exons.sort(key=lambda tup: tup[0], reverse=False)
        return coding_exons
    
    def get_introns(self):
        """ returns all introns from gene [[start, stop],*]"""
        exons = self.get_exons()
        introns = []
        last_pos = None
        for exon in exons: 
            if last_pos == None: 
                last_pos = exon[1]
            else:
                introns.append((last_pos, exon[0]))
                last_pos = exon[1]
        return introns

    def get_utr(self, which="both"): ## options are "both", "5_prime", "3_prime"
        """ returns the utrs """
        exons = self.get_exons()
        utrs = []
        utrs.append((self.txStart, self.cdsStart))
        utrs.append((self.cdsEnd, self.txEnd))
        if self.strand == "+":
            utrs = self.positive_sorted(utrs)
        elif self.strand == "-":
            utrs = self.negative_sorted(utrs)
        if which == "both":
            return self.positive_sorted(utrs)
        elif which == "3_prime":
            return utrs[1]
        elif which == "5_prime":
            return utrs[0]
        else: 
            raise Exception("which is not valid")

    def get_requested_region(self, region):
        if region == "GENIC":
            return [(self.txStart, self.txEnd)]
        elif region == "CODING":
            return self.get_coding()
        elif region == "UTR":
            return self.get_utr()
        elif region == "INTRONIC":
            return self.get_introns()    
        else:
            raise Exception("region not valid")

    ##might have to alter this 
    def get_codon_from_pos(self, pos):
        """ from a position get codon matching to that position:
        return codon, position, strand of codon of nucleotide in codon (codon, pos, strand)"""
        coding = self.get_coding()
        coding_seq = self.get_seq_from_pos(coding)
        new_pos = 0
        if self.strand == "+":
            coding = self.positive_sorted(coding)
            for exon in coding:
                if pos > exon[1]: #position is in later exon
                    new_pos = new_pos + exon[1] - exon[0] # adding length of exon
                elif exon[0] <= pos < exon[1]: #position is in this codon
                    new_pos = new_pos + pos - exon[0]
        else: 
            coding = self.negative_sorted(coding)
            for exon in coding:
                if pos < exon[0]: #position is in later exon
                    new_pos = new_pos + exon[1] - exon[0] # adding length of exon
                elif exon[0] <= pos < exon[1]: #position is in this codon
                    new_pos = new_pos + exon[1] - pos
        codon_pos = new_pos%3
        codon_start = new_pos - codon_pos 
        return (coding_seq[codon_start:codon_start+3] ,codon_pos, self.strand)

    def __str__(self):
        return "{}\t{}\t{}\t{}\t0\t{}".format(self.chrom, self.txStart, self.txEnd, self.name, self.strand)