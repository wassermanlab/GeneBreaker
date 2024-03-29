from GeneBreaker.src.variant import Variant
from GeneBreaker.src.transcript import Transcript
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
                         "Ile": ["ATT", "ATC", "ATA"], "Trp": ["TGG"], 'Met': ["ATG"]}
    aa_codes = {}
    for key, value in list(amino_acid_codons.items()):
        for code in value:
            aa_codes[code] = key
    translator = {"A": "T",
                  "T": "A",
                  "G": "C",
                  "C": "G"}

    def __init__(self, var_template, transcript: Transcript):
        """ initialize snv variant """
        try:
            Variant.__init__(self, var_template, transcript)
            self.snv_type = self.impact["SNV_TYPE"]
            self.pos = self.impact["START"]
            self.check_snv()
            if self.region == "CODING":
                if self.pos == "ANY":
                    var_dict = self.get_random_coding_SNV()
                else:
                    var_dict = self.get_directed_coding_SNV(self.pos)
            else:
                var_dict = self.get_non_coding_SNV()  # just make change to this
            if var_dict == False:
                raise Exception("Specified SNV cannot be made")
            self.pos = var_dict["pos"]
            self.ref = str(var_dict["ref"])
            self.alt = str(var_dict["alt"])
        except Exception as e:
            raise(e)

    def check_snv(self):
        if self.type != "SNV":
            raise ValueError("Must be SNV type")
        if self.snv_type not in ["MISSENSE", "NONSENSE", "SYNONYMOUS", "STOPLOSS",  "A", "T", "G", "C", "ANY"]:
            raise ValueError(
                "TYPE must be missense, nonsense, synonymous or a base")
        if self.region != "CODING" and self.snv_type in ["MISSENSE", "NONSENSE", "SYNONYMOUS", "STOPLOSS"]:
            raise ValueError("missense, nonsense, synonymous, and stoploss types are not valid for non-coding region")
        if type(self.pos) == int:
            self.pos = self.pos - 1
            self.check_location(self.pos)
        elif self.pos != "ANY":
            raise ValueError("locations must be ANY or and int")       
        if type(self.pos) == int and self.snv_type == "STOPLOSS": ## directed stoploss
            if self.transcript.get_strand() == 1:
                self.transcript.get_cdsStart()
                valid_range = range(self.transcript.get_cdsEnd()-3, self.transcript.get_cdsEnd())
            else: 
                valid_range = range(self.transcript.get_cdsStart(), self.transcript.get_cdsStart()+3)
            if self.pos not in valid_range:
                raise ValueError("directed stop loss must be at the stop codon position")

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
        missense mutation cannot be made there the alternate variant """
        alternate_codons = self.get_alternate_codons(codon, pos)
        missense_aa = self.amino_acid_codons.copy()
        missense_aa.pop(self.aa_codes[codon])
        missense_aa = [inner for outer in list(missense_aa.values())
                       for inner in outer]
        random.shuffle(alternate_codons)
        for i in alternate_codons:
            if i in missense_aa:
                return i[pos]
        return False

    def nonsense_mutation(self, codon, pos):
        """ make nonsense mutation in codon, returns False if 
        nonsense mutation cannot be made there the alternate variante """
        alternate_codons = self.get_alternate_codons(codon, pos)
        for i in alternate_codons:
            if i in self.stop_codons:
                return i[pos]
        return False

    def synonymous_mutation(self, codon, pos):
        """ make synonymous mutation in codon, returns False if 
        synonymous mutation cannot be made there the alternate variant"""
        if codon not in ["ATG", "TAA", "TAG", "TGA"]:
            alternate_codons = self.get_alternate_codons(codon, pos)
            synonymous_aa = self.amino_acid_codons[self.aa_codes[codon]]
            random.shuffle(alternate_codons)
            for i in alternate_codons:
                if i in synonymous_aa:
                    return i[pos]
            return False
        return False

    def stoploss_mutation(self, codon, pos):
        if codon not in self.stop_codons: ## if codon requested is not a stop codon
            return False
        else: 
            alternate_codons = self.get_alternate_codons(codon, pos)
            random.shuffle(alternate_codons)
            for i in alternate_codons:
                if i not in self.stop_codons:
                    return i[pos]
            return False

    def any_mutation(self, codon, pos):
        """ make random mutation in codon, returns False if 
        mutation cannot be made there the alternate variant"""
        choices = ['A', 'T', 'G', 'C']
        choices.remove(codon[pos].upper())
        alt = random.choice(choices)
        return alt

    def get_non_coding_SNV(self):
        """ make random or directed SNV """
        # get ranges
        region_range = self.get_region_range()
        if self.pos == "ANY":
            pos = random.choice(region_range)
        else:
            pos = self.pos
        shift_pos = pos - self.transcript.get_start()
        ref = self.transcript.get_seq()[shift_pos]
        if self.snv_type == "ANY":  # user has requested random base change
            choices = ['A', 'T', 'G', 'C']
            choices.remove(ref.upper())
            alt = random.choice(choices)
        else:  # user has requested specific base change
            alt = self.snv_type
        return {"pos": pos, "ref": ref, "alt": alt}

    def get_directed_coding_SNV(self, loc):
        """ get coding SNV, returns  dict{pos, ref, alt} or False"""
        codon_tuple = self.transcript.get_codon_from_pos(loc)
        codon = codon_tuple[0]
        codon_pos = codon_tuple[1]
        strand = codon_tuple[2]

        if self.snv_type == "MISSENSE":
            alt = self.missense_mutation(codon, codon_pos)
        elif self.snv_type == "NONSENSE":
            alt = self.nonsense_mutation(codon, codon_pos)
        elif self.snv_type == "SYNONYMOUS":
            alt = self.synonymous_mutation(codon, codon_pos)
        elif self.snv_type == "STOPLOSS":
            alt = self.stoploss_mutation(codon, codon_pos)
        elif self.snv_type == "ANY":
            alt = self.any_mutation(codon, codon_pos)
        else:
            alt = self.snv_type
        if alt == False:
            raise Exception("Specified " + self.snv_type + " mutation cannot be made in this position or if ANY was selected than at ANY position in this gene.")
        # strand is + or user specified mutation
        elif strand == 1:
            return {"pos": loc,
                    "ref": codon[codon_pos],
                    "alt": alt}
        elif self.snv_type in ["A", "T", "G", "C"]:
            return {"pos": loc,
                    "ref": self.translator[codon[codon_pos]],
                    "alt": alt}
        else:  # strand is -
            return {"pos": loc,
                    "ref": self.translator[codon[codon_pos]],
                    "alt": self.translator[alt]}

    def get_random_coding_SNV(self):
        """ get coding SNV, returns  dict{pos, ref, alt}"""
        # get available positions # get requested region
        if self.snv_type == "STOPLOSS": 
            if self.transcript.get_strand() == 1: 
                end = self.transcript.get_cdsEnd()
                region_range = range(end-3, end)
            else: 
                start = self.transcript.get_cdsStart() 
                region_range = range(start, start+3)
        else:
            region_range = self.get_region_range()
        # while loop randomly choosing positions in the available
        region_range = list(region_range)
        alt = False
        while len(region_range) > 0 and not alt:
            loc = random.choice(region_range)
            region_range.remove(loc)
            snv = self.get_directed_coding_SNV(loc)
            if snv != False:
                alt = True
                return snv
