from simulator.variant import Variant
from simulator.gene import Gene
import random

class SingleNucleotideVariant(Variant): 
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    amino_acid_codons = {
        "Phe": ["TTT", "TTC"],
        "Tyr": ["TAT", "TAC"],
        "His": ["CAT", "CAC"],
        "Gln": ["CAA", "CAG"],
        "Asn": ["AAT", "AAC"],
        "Lys": ["AAA", "AAG"],
        "Asp": ["GAT", "GAC"],
        "Glu": ["GAA", "GAG"],
        "Cys": ["TGT", "TGC"],
        "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "Leu": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        "Pro": ["CCT", "CCC", "CCA", "CCG"],
        "Gly": ["GGT", "GGC", "GGA", "GGG"],
        "Val": ["GTT", "GTC", "GTA", "GTG"],
        "Thr": ["ACT", "ACC", "ACA", "ACG"],
        "Ala": ["GCT", "GCC", "GCA", "GCG"],
        "Ile": ["ATT", "ATC", "ATA"],
        "Trp": ["TGG"]}
    aa_codes = {} 
    for key, value in amino_acid_codons.iteritems():
        for code in value: 
            aa_codes[code] = key

    
    def __init__(self, var_template):
        """ initialize snv variant """
        try:
            Variant.__init__(self, var_template)
            if self.type_impact not in ["MISSENSE", "NONSENSE", "SILENT", "A", "T", "G", "C", "ANY"]:
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
            region_range = region_range + range(region[0], region[1])
        if self.location == "ANY":
            pos = random.choice(region_range) 
        else: 
            if self.location not in region_range:
                raise Exception("location requested is not in region requested")
            pos = self.location
        ref = gene.get_seq()[pos]
        if self.type_impact == "ANY": # user has requested random base change
            choices = ['A', 'T', 'G', 'C'].remove(ref.upper())
            alt = random.choice(choices)
        else: # user has requested specific base change
            alt = self.type_impact
        return {"pos": pos, "ref": ref, "alt": alt}


    def get_alternate_codons(self, codon, pos):
        choices = ['A', 'T', 'G', 'C']
        choices.remove(codon[pos].upper())
        alternate_codons = []
        list_codon = list(codon.upper())
        for base in choices: 
            list_codon[pos] = base
            alternate_codons.append("".join(list_codon))
        return alternate_codons


    def missense_mutation(self, codon, pos): 
        """ make missense mutation in codon, returns False if 
        silent mutation cannot be made there the alternate variant """
        alternate_codons = self.get_alternate_codons(codon, pos)
        missense_aa = self.amino_acid_codons
        missense_aa.pop(self.aa_codes[codon])
        missense_aa = [inner for outer in missense_aa.values() for inner in outer]
        random.shuffle(alternate_codons)
        for i in alternate_codons:
            if i in missense_aa:
                return i[pos]
        return False


    def nonsense_mutation(self, codon, pos):
        """ make nonsense mutation in codon, returns False if 
        silent mutation cannot be made there the alternate variante """
        alternate_codons = self.get_alternate_codons(codon, pos)
        for i in alternate_codons:
            if i in self.stop_codons:
                return i[pos]
        return False


    def silent_mutation(self, codon, pos): 
        """ make silent mutation in codon, returns False if 
        silent mutation cannot be made there the alternate variant"""
        # Todo: code entire logic and test
        alternate_codons = self.get_alternate_codons(codon, pos)
        silent_aa = self.amino_acid_codons[self.aa_codes[codon]]
        random.shuffle(alternate_codons)
        for i in alternate_codons:
            if i in silent_aa:
                return i[pos]
        return False
    

    def any_mutation(self, codon, pos):
        """ make random mutation in codon, returns False if 
        silent mutation cannot be made there the alternate variant"""
        choices = ['A', 'T', 'G', 'C'].remove(codon[pos].upper())
        alt = random.choice(choices)
        return alt

    
    def get_codon_from_pos(self, gene, pos):
        """ from position get codon matching to position
        return codon and position of nucleotide in codon (codon, pos) """
        # Todo: code entire logic and test 
        seq = gene.get_seq()
        coding_regions = gene.get_requested_region("CODING")
        coding_regions = sorted(coding_regions, key=lambda x: x[0])
        gene = ""
        shift_pos = 0
        last_exon = None        
        coding_pos = None       # will be used to get the position requested without introns 
        for exon in coding_regions:
            gene = gene + seq[exon[0]:exon[1]]
            if last_exon is None: 
                last_exon = exon[1]
            else: 
                shift_pos = shift_pos + exon[0] - last_exon
                last_exon = exon[1]
            if exon[0] <= pos < exon[1]:
                coding_pos = pos - shift_pos
        position_in_codon = coding_pos%3
        if coding_pos < 3:
            return(seq[0:3], position_in_codon)
        else:  
            start = coding_pos - position_in_codon
            stop = start + 3
            return(seq[start:stop], position_in_codon)


    def get_random_coding_SNV(self, gene):
        """ get coding SNV, returns  dict{pos, ref, alt}"""
        # Todo: code entire logic and test
        None


    def get_directed_coding_SNV(self, gene):
        """ get coding SNV, returns  dict{pos, ref, alt}""" 
        # Todo: code entire logic and test
        None   
    

    def get_vcf_row(self, gene):
        """ get variant row tab delimitted  """
        chrom = str(gene.get_chr())
        if self.region == "CODING": 
            if self.location == "ANY":
                var_dict = self.get_random_coding_SNV(gene)
            else: 
                var_dict = self.get_directed_coding_SNV(gene)
        else: 
            var_dict = self.get_non_coding_SNV(gene)
        pos = str(var_dict["pos"])
        ref = str(var_dict["ref"])
        alt = str(var_dict["alt"])
        ID = "_".join(["snv", pos, str(self.type_impact)])
        return "\t".join([chrom, pos, ID, ref, alt])