from simulator.gene import Gene 
from simulator.variant import Variant

class Variants: 

    def __init__(self, gene_json, variants_json):
        """ 
        creates variants object which has gene, inheritance, var1, var2 
        """
        try: 
            self.gene = Gene(gene_json) # may have to expand to phenotype 
            self.inheritance = variants_json["INHERITANCE"]
            self.var1 = Variant(variants_json["VAR1"]) ## here check that it follows the general variant rules
            self.var2 = Variant(variants_json["VAR2"])
        except:
            print("Check that the variants input is correct and follows the schema")

    def variant_2_VCFR(self, variant): 
        """ turns variant template into variant, returning a dict with the 
        following keys: chrom, pos, id, ref, alt """
        ## here turn each superclass variant into its respective child class
        return {"chrom": 1, 
                "pos": 1, 
                "ID": "ABC",
                "REF": "A",
                "ALT": "C"}


    def save_vcf_output(file_name, trio)
        """ saves a vcf output of the variants
        if  trio == false then it just saves the child
        if trio == true it saves a child and 2 parents"""
        return False