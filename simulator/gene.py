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
            print("IOError: has occured, check tht gene file name exists")
        except KeyError:
            print("KeyError: check that Gene JSON is formatted correctly")
        except:
            print("Unexpected error")
        

    def get_exons(self):
        """ returns all exons from gene in the form of [[seq, start, stop],*]"""
        exons = []
        for i in self.elements:
            if i[0] == "EXON":
                exons.append([self.seq[i[1]:i[2]], i[1], i[2]])
        return exons

    def get_introns(self):
        """ returns all introns from gene [[seq, start, stop],*]"""
        introns = []
        for i in self.elements:
            if i[0] == "INTRON":
                introns.append([self.seq[i[1]:i[2]], i[1], i[2]])
        return introns

    def get_upstream_utr(self):
        """ returns the 5' utr """
        return False

    def get_downstream_utr(self):
        """ returns the 3' utr """
        return False

    def get_promoter(self):
        """ returns promoter """
        return False

    def get_enhancer(self):
        """ returns all enhancer/silencers """
        return False