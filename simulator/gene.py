from Bio.Seq import Seq
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
        return False