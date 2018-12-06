from simulator.transcript import Transcript 
from simulator.variant import Variant
from simulator.indel import Indel
from simulator.single_nucleotide_variant import SingleNucleotideVariant as SNV
from simulator.short_tandem_repeat import ShortTandemRepeat
import json

class Variants: 

    def __init__(self, variants_file):
        """ 
        creates variants object which has gene, inheritance, var1, var2 
        """
        try: 
            f = open(variants_file)
            variants_json = json.load(f)
            variants_json
            self.transcript = Transcript(variants_json["GENE_UID"])
            self.var1 = variants_json["VAR1"]
            self.var2 = variants_json["VAR2"]
            self.sex =  variants_json["SEX"]
            f.close()
        except:
            print("Check that the variants input is correct and follows the schema")

    def variants_2_VCF(self): 
        """ turns variant template into variant, string representing vcf """
        r1 = ""
        r2 = ""
        variants = [self.var1]
        if self.var2 != "NONE":
            variants.append(self.var2)
        for index, var in enumerate(variants):
            var_type = var["TYPE"]
            row = ""
            if var_type == "SNV":
                snv = SNV(var)
                row = snv.get_vcf_row(self.transcript)
            if var_type == "INDEL":
                indel = Indel(var)
                row = indel.get_vcf_row(self.transcript)
            if var_type == "STR":
                STR = ShortTandemRepeat(var)
                row = STR.get_vcf_row(self.transcript)
            if index == 0:
                r1 = row + "\n"
            if index == 1:
                r2 = row + "\n"
        return (r1, r2)


    def save_vcf_output(self, file_name):
        #Todo: impliment inheritance working
        """ saves a vcf output of the variants
        if  trio == false then it just saves the child
        if trio == true it saves a child and 2 parents"""
        vcf = self.variants_2_VCF()
        header = "##fileformat=VCFv4.2\n"
        header = header + "##fileDate=20090805\n"
        header = header + "##source=variant_simulator\n"
        header = header + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\"\n"
        header = header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND\n" 
        # TODO: implement this 
        f = open(file_name,"w+")
        f.write(header+vcf[0]+vcf[1])
        f.close()