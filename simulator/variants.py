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
            self.inheritance = variants_json["INHERITANCE"]
            self.var1 = variants_json["VAR1"]
            self.var2 = variants_json["VAR2"]
            self.trio = variants_json["TRIO"]
            f.close()
        except:
            print("Check that the variants input is correct and follows the schema")

    def variants_2_VCF(self): 
        """ turns variant template into variant, string representing vcf """
        h = "#CHROM\tPOS\tID\tREF\tALT\n"
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
        return (h, r1, r2)


    def save_vcf_output(self, file_name):
        #Todo: impliment inheritance working
        """ saves a vcf output of the variants
        if  trio == false then it just saves the child
        if trio == true it saves a child and 2 parents"""
        vcf = self.variants_2_VCF()
        if self.trio == 'SINGLE':
            f = open(file_name,"w+")
            f.write(vcf[0]+vcf[1]+vcf[2])
            f.close()
        elif self.trio == 'TRIO':
            if self.inheritance == "DE-NOVO":
                f = open(file_name+".child","w+")
                f.write(vcf[0]+vcf[1]+vcf[2])
                f.close()
                f = open(file_name+".mother","w+")
                f.write(vcf[0])
                f.close()
                f = open(file_name+".father","w+")
                f.write(vcf[0])
                f.close()
            elif self.inheritance == "BI-PARENTAL":
                if vcf[2] == "":
                    f = open(file_name+".child","w+")
                    f.write(vcf[0]+vcf[1])
                    f.close()
                    f = open(file_name+".mother","w+")
                    f.write(vcf[0]+vcf[1])
                    f.close()
                    f = open(file_name+".father","w+")
                    f.write(vcf[0]+vcf[1])
                    f.close()
                else:
                    f = open(file_name+".child","w+")
                    f.write(vcf[0]+vcf[1]+vcf[2])
                    f.close()
                    f = open(file_name+".mother","w+")
                    f.write(vcf[0]+vcf[1])
                    f.close()
                    f = open(file_name+".father","w+")
                    f.write(vcf[0]+vcf[2])
                    f.close()
            elif self.inheritance == "MATERNAL":
                f = open(file_name+".child","w+")
                f.write(vcf[0]+vcf[1]+vcf[2])
                f.close()
                f = open(file_name+".mother","w+")
                f.write(vcf[0]+vcf[1])
                f.close()
                f = open(file_name+".father","w+")
                f.write(vcf[0])
                f.close()
            elif self.inheritance == "PATERNAL":
                f = open(file_name+".child","w+")
                f.write(vcf[0]+vcf[1]+vcf[2])
                f.close()
                f = open(file_name+".mother","w+")
                f.write(vcf[0])
                f.close()
                f = open(file_name+".father","w+")
                f.write(vcf[0]+vcf[1])
                f.close()
