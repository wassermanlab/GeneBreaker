from Bio.Seq import Seq


class Gene:

    def __init__(self, seq, chr, start_pos, stop_pos, elements):
        """ create gene object with sequence, chromosome, start position, 
        stop position and elements. 
        elements tuples in this form (type, start, stop) where type can be one 
        of the following: exon, intron, utr, promoter, enhancer
        """
        self.seq = seq
        self.chr = chr
        self.start_reg = start_pos
        self.stop_reg = stop_pos
        self.elements = []
        for el in elements:
            self.elements.append(el)

    def get_exons():
        """ returns all exons from gene """
        return False

    def get_introns():
        """ returns all introns from gene """
        return False

    def get_upstream_utr():
        """ returns the 5' utr """
        return False

    def get_downstream_utr():
        """ returns the 3' utr """
        return False

    def get_promoter():
        """ returns promoter """
        return False

    def get_enhancer():
        """ returns all enhancer/silencers """
        