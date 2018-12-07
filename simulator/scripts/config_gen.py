# todo
# fix TRIO responses
import argparse
import json
from GUD2.ORM import Gene
from GUD2.ORM import ShortTandemRepeat
from GUD2.ORM import ClinVar
import os
import sys
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
import warnings
from simulator.transcript import Transcript

# Establish a SQLalchemy session w/ GUD
db_name = "mysql://{}:@{}:{}/{}".format("ontarget_r", ## change to hg19
                                        "ontarget.cmmt.ubc.ca", "5506", "tamar_test") 
try:
    engine = create_engine(db_name, echo=False)
    session = Session(engine) 
except:
    raise ValueError("Cannot connect to GUD: %s" % "hg19")

def get_str_impact(transcript, var_region):
    regions = transcript.get_requested_region(var_region)
    STR = ShortTandemRepeat()
    chrom = transcript.chrom
    count = 0 
    str_list = []

    for region in regions:
        str_list = str_list + STR.select_by_location(session, chrom, region[0], region[1])
    
    if len(str_list) == 0:
        print "There is no STR in that region please select another variant"
        return False
    print "UID\tMOTIF\tPATHOGENICITY"
    for s in str_list:
        print str(s[0].uid) + "\t" + str(s[0].motif) + "\t" + str(s[0].pathogenicity) 

    uid = int(raw_input("""Input the uid of your STR: """))
    STR = STR.select_by_uid(session, uid)
    str_length = int(raw_input("""input the size of your short tandem repeat, negative numbers are deletions positive numbers are insersions: """))
    var_impact = {"CHROM": STR[1].chrom, "START": STR[1].start, "END": STR[1].end, "STR": str_length}
    return var_impact

def get_snv_impact(var_region):
    if var_region == "CODING":
        var_impact = str(raw_input("""Input the impact of your SNV, valid inputs are 'A', 'T', 'G', 'C', 'ANY', 'MISSENSE', 'NONSENSE' or 'SYNONYMOUS': """))
        if var_impact not in ['A', 'T', 'G', 'C', 'ANY', 'MISSENSE', 'NONSENSE', 'SYNONYMOUS']:
            raise Exception("Not valid impact type for CODING SNV")
    else:
        var_impact = str(raw_input("""Input the impact of your SNV, valid inputs are 'A', 'T', 'G', 'C', or 'ANY': """))
        if var_impact not in ['A', 'T', 'G', 'C', 'ANY']:
            raise Exception("Not valid impact type for SNV")
    return var_impact

def get_indel_impact():
    var_impact = int(raw_input("""input the size of your indel, negative numbers are deletions positive numbers are insersions: """))
    return var_impact

def get_clinvar_impact(transcript, var_region):
    print "HERE"
    regions = transcript.get_requested_region(var_region)
    clinvar = ClinVar()
    chrom = transcript.chrom
    count = 0 
    clinvar_list = []
    print regions
    for region in regions:
        clinvar_list = clinvar_list + clinvar.select_by_location(session, chrom, region[0], region[1])
    
    print clinvar_list
    if len(clinvar_list) == 0:
        print "There are no ClinVar variants in that region please select another variant"
        return False
    print "clinvarID\tpos(1-based)\tref\talt\tCLNSIG\tCLNDN"
    for s in str_list:
        clinvar = s[0]
        region = s[1]
        print '\t'.join([str(clinvar.clinvarID), str(region.start+1), str(clinvar.ref), 
        str(clinvar.alt), str(clinvar.CLNSIG), str(clinvar.CLNDN)])

    uid = int(raw_input("""Input the uid of your STR: """))
    clinvar = clinvar.select_by_uid(session, uid)
    var_impact = clinvar.clinvarID
    return var_impact


def get_variant_dict(transcript, sex, variant = "var1"):  
    if transcript.chrom in ["chrX", "chrY"] and sex == "XY-MALE":
        zygosity_options = ["HEMIZYGOUS"]
    elif variant == "var2":
        zygosity_options = ["HETEROZYGOUS"]  
    else:
        zygosity_options = ["HOMOZYGOUS", "HETEROZYGOUS"]  
    
    var_impact = False
    while var_impact == False:
        var_dict = dict()
        var_type = str(raw_input("""What type of variant would you like, options are 'SNV', 'INDEL', 'STR', 'ClinVar': """))  # add other variant types
        var_region = str(raw_input("""In what region would you like your variant to be, options are CODING', 'UTR', 'INTRONIC': """))  # todo add promoter and enhancer
        var_zygosity = str(raw_input("What zygosity would you like your variant to have, options are " + str(zygosity_options) + ": "))
        var_location = None
        if var_type == "INDEL":
            var_impact = get_indel_impact()
        elif var_type == "SNV":
            var_impact = get_snv_impact(var_region)
        elif var_type == "STR":
            var_impact = get_str_impact(transcript, var_region)       
            var_location = "NONE"
        elif var_type == "ClinVar":
            var_impact = get_clinvar_impact(transcript, var_region)
            var_location = "NONE"
        else:
            raise Exception("variant type specified is not valid")
    
    if var_location != "NONE": 
        var_location = str(
            raw_input("Specify a 1 based location or input 'ANY': "))
        if var_location.isdigit():
            var_location = long(var_location) - 1

    if transcript.chrom in ["chrX", "chrY"] and sex == "XY-MALE" and var_zygosity != "HEMIZYGOUS":
        raise Exception("not valid zygosity.")
    elif var_zygosity not in ["HOMOZYGOUS", "HETEROZYGOUS", "HEMIZYGOUS"]:
        raise Exception("not valid zygosity.")
    var_dict["TYPE"] = var_type
    var_dict["REGION"] = var_region
    var_dict["IMPACT"] = var_impact
    var_dict["LOCATION"] = var_location
    var_dict["ZYGOSITY"] = var_zygosity
    return var_dict

def get_main_dict():
    config_dict = dict()
    # highly coupled with GUD, need to figure out a better way to do this
    
    gene_sym = str(
        raw_input("Enter the gene symbol of the gene you would like to use: "))
    gene = Gene()
    genes = gene.select_by_name(session, gene_sym)  # GALK2
    if len(genes) == 0:
        raise Exception("No genes by that gene symbol")
    print "Gene_UID\ttranscript_start\ttranscript_end"
    for g in genes:
        full_entry = gene.select_by_uid_joined(session, g.uid)
        print str(full_entry[0].uid) + "\t" + str(full_entry[1].start) + "\t" + str(full_entry[1].end)
    gene_uid = int(
        raw_input("what is the UID of the transcript would you like to use: "))
    gene = gene.select_by_uid_joined(session, gene_uid)  # get gene that we need
    config_dict["GENE_NAME"] = gene[0].name2
    config_dict["GENE_UID"] = gene_uid
    
    sex = str(raw_input("What sex would you like the proband, options are 'XX-FEMALE' or 'XY-MALE': "))
    config_dict["SEX"] = sex
    if sex not in ['XX-FEMALE', 'XY-MALE']:
        raise Exception("Not valid sex")

    return config_dict

def main():
    config_name = str(raw_input("enter the name of your config file: "))
    config_dict = get_main_dict()
    gene_uid = config_dict["GENE_UID"]
    sex = config_dict["SEX"]
    transcript = Transcript(gene_uid)
    if sex == "XX-FEMALE" and transcript.chrom == "chrY":
        raise Exception("Cannot have a variant on the y chromosome in a female.")
    # get variant I
    print("""===========Variant I===========""")
    config_dict["VAR1"] = get_variant_dict(transcript, sex)
    
    # get variant II
    if config_dict["VAR1"]["ZYGOSITY"] == "HETEROZYGOUS":
        print("""===========Variant II===========""")
        var2_ans = str(raw_input("Would you like a second variant [y/n]: "))
        if var2_ans not in ["y", "n"]:
            raise Exception("Not valid response.")
        elif var2_ans == "y":
            config_dict["VAR2"] = get_variant_dict(transcript, sex, "var2")
        else:
            config_dict["VAR2"] = "NONE"
    else: 
        config_dict["VAR2"] = "NONE"
    
    # write the file
    if config_name.endswith(".json"):
        with open(config_name, 'w+') as fp:
            json.dump(config_dict, fp, sort_keys=True, indent=4)
    else:
        with open(config_name + ".json", 'w+') as fp:
            json.dump(config_dict, fp, sort_keys=True, indent=4)

if __name__ == "__main__":
    main()
