from simulator.transcript import Transcript 
from simulator.variant import Variant
from simulator.indel import Indel
from simulator.single_nucleotide_variant import SingleNucleotideVariant as SNV
from simulator.short_tandem_repeat import ShortTandemRepeat
from simulator.clinvar import ClinVar
from simulator.mei import MEI
from simulator.copy_number_variant import CopyNumberVariant
import json
from datetime import date

class Variants: 

    def __init__(self, variants_file):
        """ 
        creates variants object which has gene, inheritance, var1, var2 
        """ 
        f = open(variants_file)
        variants_json = json.load(f)
        variants_json
        self.transcript = Transcript(variants_json["GENE_UID"])
        self.var1 = variants_json["VAR1"]
        self.var2 = variants_json["VAR2"]
        self.sex =  variants_json["SEX"]
        f.close()

    def variants_2_VCF(self, format="simple"): 
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
            if var_type == "MEI":
                mei = MEI(var)
                row = mei.get_vcf_row(self.transcript)
            if var_type == "ClinVar":
                cln = ClinVar(var)
                row = cln.get_vcf_row(self.transcript)
            if var_type == "CNV":
                cnv = CopyNumberVariant(var)
                if format == "simple":
                    row = cnv.get_vcf_row(self.transcript, "simple")
                if format == "4.2":
                    row = cnv.get_vcf_row(self.transcript, "4.2")
            if index == 0:
                r1 = row + "\n"
            if index == 1:
                r2 = row + "\n"
        return (r1, r2)


    def save_vcf_output(self, file_name):
        variant_types = [self.var1["TYPE"]]
        if self.var2 != "NONE":
            variants.append(self.var2[["TYPE"]])
        
        if "CNV" in variant_types:
            ## make complex VCF
            vcf = self.variants_2_VCF(format="4.2")
            header = "##fileformat=VCFv4.2\n"
            header = header + "##fileDate=" + str(date.today()) + "\n"
            header = header + "##source=variant_simulator\n"
            header = header + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\"\n"
            header = header + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\"\n"
            header = header + "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\"\n"
            header = header + "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\"\n"
            header = header + "##ALT=<ID=DUP,Description=\"Duplication\"\n"
            header = header + "##ALT=<ID=DEL,Description=\"Deletion\"\n"
            header = header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND\n" 
            # TODO: implement this 
            f = open("vcf4.2."+file_name,"w+")
            f.write(header+vcf[0]+vcf[1])
            f.close()
    
        vcf = self.variants_2_VCF(format="simple")
        header = "##fileDate=" + str(date.today()) + "\n"
        header = header + "##source=variant_simulator\n"
        header = header + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\"\n"
        header = header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND\n" 
        # TODO: implement this 
        f = open(file_name,"w+")
        f.write(header+vcf[0]+vcf[1])
        f.close()