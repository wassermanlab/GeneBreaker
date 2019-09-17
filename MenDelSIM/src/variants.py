from MenDelSIM.src.transcript import Transcript
from MenDelSIM.src.variant import Variant
from MenDelSIM.src.clinvar import ClinVar
from MenDelSIM.src.copy_number_variant import CopyNumberVariant
from MenDelSIM.src.indel import Indel
from MenDelSIM.src.mei import MEI
from MenDelSIM.src.single_nucleotide_variant import SingleNucleotideVariant as SNV
from MenDelSIM.src.short_tandem_repeat import ShortTandemRepeat
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

    def save_vcf_output(self):
        #     # TODO: implement this
        # if "CNV" in variant_types:
        #     # make complex VCF
        #     vcf = self.variants_2_VCF(format="4.2")
        #     header = "##fileformat=VCFv4.2\n"
        #     header = header + "##fileDate=" + str(date.today()) + "\n"
        #     header = header + "##source=variant_simulator\n"
        #     header = header + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\"\n"
        #     header = header + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\"\n"
        #     header = header + "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\"\n"
        #     header = header + "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\"\n"
        #     header = header + "##ALT=<ID=DUP,Description=\"Duplication\"\n"
        #     header = header + "##ALT=<ID=DEL,Description=\"Deletion\"\n"
        #     header = header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND\n"
        #     f = open("vcf4.2."+file_name, "w+")
        #     f.write(header+vcf[0]+vcf[1])
        #     f.close()

        vcf = self.variants_2_VCF()
        return {"var1": vcf[0],
                "var2": vcf[1]}
