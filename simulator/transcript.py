import json
import argparse
import json
from GUD.ORM import Gene
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
from lxml import etree
from Bio.Seq import Seq

class Transcript:

    def __init__(self, name, index):
        """ create transcript object 
        """
        # Establish a SQLalchemy session w/ GUD
        db_name = "mysql://{}:@{}:{}/{}".format("ontarget_r",
            "ontarget.cmmt.ubc.ca", "5506", "hg19")

        try:
            self.name = name
            engine = create_engine(db_name, echo=False)
            session = Session(engine)
            gene = Gene()
            transcripts = Gene.select_by_name(session, self.name)
            gene = transcripts[index]
            self.cdsStart   = long(gene.cdsStart)
            self.cdsEnd     = long(gene.cdsEnd)
            self.exonStarts = gene.exonStarts
            self.exonEnds   = gene.exonEnds
            self.chrom      = gene.chrom
            self.txStart    = long(gene.txStart)
            self.txEnd      = long(gene.txEnd)
            self.strand     = gene.strand
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
            return Seq(cut).reverse_complement()

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
        """ returns all coding regions from gene in the form of [[start, stop],*]"""
        exons = self.get_exons()
        coding_exons = []
        if self.cdsStart == self.cdsEnd == self.txEnd:
            return coding_exons
        for exon in exons: 
            if exon[0] >= self.cdsStart and exon[1] <= self.cdsEnd: # exon is fully within coding region
                coding_exons.append(exon)
            elif exon[0] < self.cdsStart and exon[1] <= self.cdsEnd: #exon start is out of range but exon end is withing range
                coding_exons.append((self.cdsStart, exon[1]))
            elif exon[0] >= self.cdsStart and exon[1] > self.cdsEnd: #exon end is out of range but start is within range
                coding_exons.append((exon[0], self.cdsEnd))
            elif exon[0] < self.cdsStart and exon[1] > self.cdsEnd: #exon end and start are out of range 
                coding_exons.append((self.cdsStart, self.cdsEnd))
            else: #exon is completely out of range
                None 
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
        if region == "CODING":
            return self.get_coding()
        elif region == "UTR":
            return self.get_utr()
        elif region == "INTRONIC":
            return self.get_introns()    
        else:
            raise Exception("region not valid")

    ##might have to alter this 
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