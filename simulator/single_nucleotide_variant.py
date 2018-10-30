from simulator.variant import Variant
from simulator.gene import Gene
import random

class SingleNucleotideVariant(variant.Variant): 
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    amino_acid_codons = [("Phe", "TTT"), ("Phe", "TTC"),
    ("Leu", "TTA"), ("Leu","TTG"), ("Leu","CTT"), ("Leu","CTC"), ("Leu","CTA"), ("Leu","CTG"),
    ("Val", "GTT"), ("Val", "GTC"), ("Val", "GTA"), ("Val", "GTG"),
    ("Ser", "TCT"), ("Ser", "TCC"), ("Ser", "TCA"), ("Ser", "TCG"), ("Ser", "AGT"), ("Ser", "AGC"),
    ("Pro", "CCT"), ("Pro", "CCC"), ("Pro", "CCA"), ("Pro", "CCG"),
    ("Thr", "ACT"), ("Thr", "ACC"), ("Thr", "ACA"), ("Thr", "ACG"),
    ("Ala", "GCT"), ("Ala", "GCC"), ("Ala", "GCA"), ("Ala", "GCG"),
    ("Tyr", "TAT"), ("Tyr", "TAC"),
    ("His", "CAT"), ("His", "CAC"),
    ("Gln", "CAA"), ("Gln", "CAG"),
    ("Asn", "AAT"), ("Asn", "AAC"),
    ("Lys", "AAA"), ("Lys", "AAG"),
    ("Asp", "GAT"), ("Asp", "GAC"),
    ("Glu", "GAA"), ("Glu", "GAG"),
    ("Cys", "TGT"), ("Cys", "TGC"),
    ("Trp", "TGG"),
    ("Arg", "CGT"), ("Arg", "CGC"), ("Arg", "CGA"), ("Arg", "CGG"), ("Arg", "AGA"), ("Arg", "AGG").
    ("Gly", "GGT"), ("Gly", "GGC"), ("Gly", "GGA"), ("Gly", "GGG").
    ("Ile", "ATT"), ("Ile", "ATC"), ("Ile", "ATA").]

    # {"TYPE": "SNV",
    #     "REGION": "CODING",
    #     "IMPACT":{"TYPE_IMPACT": MISSENSE | NONSENSE | SILENT | ATGC, "LOCATION": "ANY" | NUMBER}}

    def __init__(self, var_template):
        """ initialize snv variant """
        try:
            Variant.__init__(self, var_template)
            if self.type_impact not in ["MISSENSE", "NONSENSE", "SILENT", "A", "T", "G", "C"]:
                raise Exception("""TYPE must be missense, nonsense, silent of a base""")
            if self.type != "SNV":
                raise Exception("Must be indel type")
        except:
            print('check that type is indel and that impact is an int')
    
    def get_non_coding_SNV(self, gene):
        """ make random or directed SNV """
        # get requested region
        regions = gene.get_requested_region(self.region)
        if len(regions) == 0:
            raise Exception("region requested does not exists")
        # get ranges
        region_range = []
        for region in regions: # range must be cut so that there is no overlap
        region_range = region_range + \
        range(region[0], region[1])
        if self.location == "ANY":
            pos = random.choice(region_range) 
        else: 
            if self.location not in region_range:
                raise Exception("location requested is not in region requested")
            pos = self.location
        ref = gene.get_seq()[pos]
        choices = ['A', 'T', 'G', 'C'].remove(ref.upper())
        alt = random.choice(choices)
        return {"pos": pos, "ref": ref, "alt": alt}

    def missense_mutation(self): 
        """ make a missense mutation in random or directed place """
    
    def nonsense_mutation(self):
        """ make nonsense in random or directed place """
    
    def silent_mutation(self): 
        """ make silent mutation in random or directed place """
    
    def get_coding_SNV(self, gene):
        """ get coding SNV """
        if self.location == "ANY":
            
        else: 
    
    def get_vcf_row(self, gene):
        """ get variant row tab delimitted  """
        chrom = str(gene.get_chr())
        if self.REGION == "CODING": 
            var_dict = self.get_coding_SNV(gene)
        else: 
            var_dict = self.get_non_coding_SNV(gene)
        pos = str(var_dict["pos"])
        ref = str(var_dict["ref"])
        alt = str(var_dict["alt"])
        ID = "_".join(["snv", pos, str(self.type_impact)])
        return "\t".join([chrom, pos, ID, ref, alt])