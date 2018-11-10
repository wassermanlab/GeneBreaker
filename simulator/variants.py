from simulator.transcript import Transcript 
from simulator.variant import Variant
from simulator.indel import Indel
from simulator.single_nucleotide_variant import SingleNucleotideVariant as SNV
import json

class Variants: 

    def __init__(self, variants_file):
        """ 
        creates variants object which has gene, inheritance, var1, var2 
        """
        try: 
            f = open(variants_file)
            variants_json = json.load(f)
            self.transcript = Transcript(variants_json["GENE_NAME"], variants_json["GENE_INDEX"])
            self.inheritance = variants_json["INHERITANCE"]
            self.var1 = variants_json["VAR1"]
            self.var2 = variants_json["VAR2"]
            self.trio = variants_json["TRIO"]
            f.close()
        except:
            print("Check that the variants input is correct and follows the schema")

    def variants_2_VCF(self): 
        """ turns variant template into variant, string representing vcf """
        vcf = "#CHROM\tPOS\tID\tREF\tALT\n"
        variants = [self.var1]
        if self.var2 != "NONE":
            variants.append(self.var2)
        for var in variants:
            var_type = var["TYPE"]
            row = ""
            if var_type == "SNV":
                snv = SNV(var)
                row = snv.get_vcf_row(self.transcript)
            if var_type == "INDEL":
                indel = Indel(var)
                row = indel.get_vcf_row(self.transcript)
            vcf = vcf + row + "\n"
        return vcf


    def save_vcf_output(self, file_name):
        #Todo: impliment inheritance working
        """ saves a vcf output of the variants
        if  trio == false then it just saves the child
        if trio == true it saves a child and 2 parents"""
        vcf = self.variants_2_VCF()
        f = open(file_name,"w+")
        f.write(vcf)
        f.close()