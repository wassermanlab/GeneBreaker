from simulator.variant import Variant
from simulator.gene import Gene 

class Indel(Variant): 
    # assume var_template is of type dict already
    def __init__(self, var_template):
        try: 
            Variant.__init__(self, var_template)
            if type(self.impact) is not int:
                raise Exception("""For indel type the program expects an int 
                representing the amount to delete or insert, did not get that input.""")
            if self.type is not "INDEL":
                raise Exception("Must be indel type")
        except:
            print('check that type is indel and that impact is an int')
    

    def get_vcf_row(self, gene):
        return False
