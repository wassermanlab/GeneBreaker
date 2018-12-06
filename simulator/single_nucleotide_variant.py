from simulator.variant import Variant
from simulator.transcript import Transcript
import random

class SingleNucleotideVariant(Variant):
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    amino_acid_codons = {"Phe": ["TTT", "TTC"], "Tyr": ["TAT", "TAC"],
                         "His": ["CAT", "CAC"], "Gln": ["CAA", "CAG"], "Asn": ["AAT", "AAC"],
                         "Lys": ["AAA", "AAG"], "Asp": ["GAT", "GAC"], "Glu": ["GAA", "GAG"],
                         "Cys": ["TGT", "TGC"], "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
                         "Leu": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
                         "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
                         "Pro": ["CCT", "CCC", "CCA", "CCG"],
                         "Gly": ["GGT", "GGC", "GGA", "GGG"],
                         "Val": ["GTT", "GTC", "GTA", "GTG"],
                         "Thr": ["ACT", "ACC", "ACA", "ACG"],
                         "Ala": ["GCT", "GCC", "GCA", "GCG"],
                         "Ile": ["ATT", "ATC", "ATA"], "Trp": ["TGG"]}
    aa_codes = {}
    for key, value in amino_acid_codons.iteritems():
        for code in value:
            aa_codes[code] = key
    translator = {"A": "T",
                  "T": "A",
                  "G": "C",
                  "C": "G"}

    def __init__(self, var_template):
        """ initialize snv variant """
        try:
            Variant.__init__(self, var_template)
            if self.impact not in ["MISSENSE", "NONSENSE", "SILENT", "A", "T", "G", "C", "ANY"]:
                raise Exception(
                    """TYPE must be missense, nonsense, silent of a base""")
            if self.type != "SNV":
                raise Exception("Must be SNV type")
            if self.region != "CODING" and self.impact in ["MISSENSE", "NONSENSE", "SILENT"]:
                raise Exception("type impact not valid for non coding region")
        except:
            print('''check that type is SNV and that impact is a one of: 
            "MISSENSE", "NONSENSE", "SILENT", "A", "T", "G", "C", "ANY"''')

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
        if codon not in ["ATG", "TAA", "TAG", "TGA"]:
            alternate_codons = self.get_alternate_codons(codon, pos)
            missense_aa = self.amino_acid_codons.copy()
            missense_aa.pop(self.aa_codes[codon])
            missense_aa = [inner for outer in missense_aa.values()
                           for inner in outer]
            random.shuffle(alternate_codons)
            for i in alternate_codons:
                if i in missense_aa:
                    return i[pos]
            return False
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
        if codon not in ["ATG", "TAA", "TAG", "TGA"]:
            alternate_codons = self.get_alternate_codons(codon, pos)
            silent_aa = self.amino_acid_codons[self.aa_codes[codon]]
            random.shuffle(alternate_codons)
            for i in alternate_codons:
                if i in silent_aa:
                    return i[pos]
            return False
        return False

    def any_mutation(self, codon, pos):
        """ make random mutation in codon, returns False if 
        silent mutation cannot be made there the alternate variant"""
        choices = ['A', 'T', 'G', 'C']
        choices.remove(codon[pos].upper())
        alt = random.choice(choices)
        return alt

    def get_non_coding_SNV(self, transcript):
        """ make random or directed SNV """
        # get requested region
        regions = transcript.get_requested_region(self.region)
        if len(regions) == 0:
            raise Exception("region requested does not exists")
        # get ranges
        region_range = []
        for region in regions:  # range must be cut so that there is no overlap
            region_range = region_range + range(region[0], region[1])
        if self.location == "ANY":
            pos = random.choice(region_range)
        else:
            if self.location not in region_range:
                raise Exception(
                    "location requested is not in region requested")
            pos = self.location
        shift_pos = pos - transcript.get_start()
        ref = transcript.get_seq()[shift_pos]
        if self.impact == "ANY":  # user has requested random base change
            choices = ['A', 'T', 'G', 'C'].remove(ref.upper())
            alt = random.choice(choices)
        else:  # user has requested specific base change
            alt = self.impact
        return {"pos": pos, "ref": ref, "alt": alt}

    def get_directed_coding_SNV(self, transcript, loc):
        """ get coding SNV, returns  dict{pos, ref, alt} or False"""
        codon_tuple = transcript.get_codon_from_pos(loc)
        codon = codon_tuple[0]
        codon_pos = codon_tuple[1]
        strand = codon_tuple[2]
        if self.impact == "MISSENSE":
            alt = self.missense_mutation(codon, codon_pos)
        elif self.impact == "NONSENSE":
            alt = self.nonsense_mutation(codon, codon_pos)
        elif self.impact == "SILENT":
            alt = self.silent_mutation(codon, codon_pos)
        elif self.impact == "ANY":
            alt = self.any_mutation(codon, codon_pos)
        else:
            alt = self.impact
        if alt is False:
            return False
        # strand is + or user specified mutation
        elif strand == "+":
            return {"pos": loc,
                    "ref": codon[codon_pos],
                    "alt": alt}
        elif self.impact in ["A", "T", "G", "C"]:
            return {"pos": loc,
                    "ref": self.translator[codon[codon_pos]],
                    "alt": alt}
        else:  # strand is -
            return {"pos": loc,
                    "ref": self.translator[codon[codon_pos]],
                    "alt": self.translator[alt]}

    def get_random_coding_SNV(self, transcript):
        """ get coding SNV, returns  dict{pos, ref, alt}"""
        # get available positions # get requested region
        regions = transcript.get_requested_region(self.region)
        if len(regions) == 0:
            raise Exception("region requested does not exists")
        region_range = []
        for region in regions:  # range must be cut so that there is no overlap
            region_range = region_range + range(region[0], region[1])
        # while loop randomly choosing positions in the available
        alt = False
        while len(region_range) > 0:
            loc = random.choice(region_range)
            region_range.remove(loc)
            snv = self.get_directed_coding_SNV(transcript, loc)
            if snv is not False:
                alt = True
                return snv
                break

    def get_vcf_row(self, transcript):
        """ get variant row tab delimitted  """
        chrom = str(transcript.get_chr())
        if self.region == "CODING":
            if self.location == "ANY":
                var_dict = self.get_random_coding_SNV(transcript)
            else:
                var_dict = self.get_directed_coding_SNV(
                    transcript, self.location)
        else:
            var_dict = self.get_non_coding_SNV(transcript)
        if var_dict is False:
            raise Exception("Specified SNV cannot be made")
        pos = str(var_dict["pos"] + 1)  # changing to 1 based
        ref = str(var_dict["ref"])
        alt = str(var_dict["alt"])
        ID = "_".join(["snv", pos, str(self.impact)])
        if self.zygosity == "HOMOZYGOUS":
            zygosity = "1/1"  
        if self.zygosity == "HEMIZYGOUS":
            zygosity = "1"
        if self.zygosity == "HETEROZYGOUS":
            zygosity = "0/1"  
        return "\t".join([chrom, pos, ID, ref, alt, ".", ".", ".", "GT", zygosity])
