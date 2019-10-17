from GeneBreaker.src.transcript import Transcript
from GeneBreaker.src.variant import Variant
from GeneBreaker.src.clinvar import ClinVar
from GeneBreaker.src.copy_number_variant import CopyNumberVariant
from GeneBreaker.src.indel import Indel
from GeneBreaker.src.mei import MEI
from GeneBreaker.src.single_nucleotide_variant import SingleNucleotideVariant as SNV
from GeneBreaker.src.short_tandem_repeat import ShortTandemRepeat
import json
from datetime import date


class Variants:

    def __init__(self, variants_file):
        """ 
        creates variants object which has gene, inheritance, var1, var2 
        """
        if type(variants_file) != dict:
            f = open(variants_file)
            variants_json = json.load(f)
            f.close()
        else:
            variants_json = variants_file

        try:
            self.transcript = Transcript(
                variants_json["GENE_UID"], variants_json["GENOME"])
            self.var1 = self.make_variant(
                variants_json["VAR1"], self.transcript)
            self.var2 = self.make_variant(
                variants_json["VAR2"], self.transcript)
            self.sex = variants_json["SEX"]
            self.check_sex_zygosity()
        except Exception as e:
            raise Exception(e)

    def make_variant(self, var, transcript):
        if var == "None":
            return None
        var_type = var["TYPE"]
        if var_type == "SNV":
            variant = SNV(var, transcript)
        elif var_type == "INDEL":
            variant = Indel(var, transcript)
        elif var_type == "STR":
            variant = ShortTandemRepeat(var, transcript)
        elif var_type == "MEI":
            variant = MEI(var, transcript)
        elif var_type == "CLINVAR":
            variant = ClinVar(var, transcript)
        elif var_type == "CNV":
            variant = CopyNumberVariant(var, transcript)
        elif var_type == "CLINGEN":
            variant = CopyNumberVariant(var, transcript)
        else:
            variant = Variant(var, transcript)
        return variant

    def check_sex_zygosity(self):
        if self.sex not in ['XX', 'XY']:
            raise ValueError('SEX must be one of: XX, XY')
        if self.transcript.chrom == 'chrY' and self.sex != 'XY':
            raise ValueError('SEX of proband must be XY when selecting Y gene')
        if self.var1.zygosity in ['HOMOZYGOUS', 'HEMIZYGOUS'] and self.var2 != None:
            raise ValueError(
                'Cannot have a second variant if the first is HOMOZYGOUS or HEMIZYGOUS')
        if self.var1.zygosity != 'HEMIZYGOUS' and self.sex == 'XY' and self.transcript.get_chr() in ["chrX", "chrY"]:
            raise ValueError('With XY sex zygosity must be HEMIZYGOUS')

    def variants_2_VCF(self):
        """ turns variant template into variant, string representing vcf """
        # print(self.var2)
        r2 = ""
        if self.var2 != None:
            r2 = self.var2.get_vcf_row() 
        r1 = self.var1.get_vcf_row() 
        return (r1, r2)

    def get_variant_rows(self):
        vcf = self.variants_2_VCF()
        return {"var1": vcf[0],
                "var2": vcf[1]}
