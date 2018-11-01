import json

class Gene:

    def __init__(self, gene_json):
        """ create gene object with sequence, chromosome, start position, 
        stop position and elements. 
        elements tuples in this form (type, start, stop) where type can be one 
        of the following: exon, intron, utr, promoter, enhancer
        """
        try:
            f = open(gene_json)
            data = json.load(f)
            self.seq = data["SEQ"]
            self.chr = data["CHR"]
            self.start_reg = data["START_REG"]
            self.stop_reg = data["STOP_REG"]
            self.elements = data["ELEMENTS"]
            f.close()
        except IOError:
            print("IOError: has occured, check that gene file name exists")
        except KeyError:
            print("KeyError: check that Gene JSON is formatted correctly")
        except:
            print("Unexpected error")

    def get_seq(self): 
        """ returns all sequence from gene""" 
        return self.seq

    def get_chr(self): 
        """ returns chr from gene""" 
        return self.chr       

    def get_exons(self):
        """ returns all exons from gene in the form of [[seq, start, stop],*]"""
        exons = []
        for i in self.elements:
            if i[0] == "EXON":
                exons.append([i[1], i[2]])
        return exons

    def get_introns(self):
        """ returns all introns from gene [[seq, start, stop],*]"""
        introns = []
        for i in self.elements:
            if i[0] == "INTRON":
                introns.append([i[1], i[2]])
        return introns

    def get_utr(self):
        """ returns the utrs """
        return False

    def get_promoter(self):
        """ returns promoter """
        return False

    def get_enhancer(self):
        """ returns all enhancer/silencers """
        return False

    def get_requested_region(self, region):
        if region == "CODING":
            return self.get_exons()
        elif region == "UTR":
            return self.get_utr()
        elif region == "INTRONIC":
            return self.get_introns()
        elif region == "PROMOTER":
            return self.get_promoter()
        elif region == "ENHANCER":
            return self.get_enhancer()    
        else:
            raise Exception("region not valid")

    def get_codon_from_pos(self, pos):
        """ from position get codon matching to position
        return codon and position of nucleotide in codon (codon, pos) """
        # Todo: code entire logic and test
        coding_regions = self.get_exons()
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
            return(self.seq[start:stop], position_in_codon)