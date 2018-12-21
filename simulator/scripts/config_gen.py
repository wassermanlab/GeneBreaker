# todo
# fix TRIO responses
import argparse
import json
import re
import os
import sys
from GUD2.ORM import Gene
from GUD2.ORM import ShortTandemRepeat
from GUD2.ORM import ClinVar
from GUD2.ORM import Chrom
from GUD2.ORM import CNV
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
from simulator.transcript import Transcript

# Establish a SQLalchemy session w/ GUD
db_name = "mysql://{}:@{}:{}/{}".format("ontarget_r", ## change to hg19
                                        "ontarget.cmmt.ubc.ca", "5506", "tamar_test") 
try:
    engine = create_engine(db_name, echo=False)
    session = Session(engine) 
except:
    raise ValueError("Cannot connect to GUD: %s" % "hg19")

def parse_region(transcript, var_region):
    """returns a tuple (chrom, (start, stop))"""
    if re.match("^chr([XY]|[1-9]|1[0-9]|2[0-2]):\d+-\d+$", var_region) is None:
        raise Exception("custom region is not in correct format")
    chrom = var_region.split(":")[0]
    start = int(var_region.split(":")[1].split("-")[0]) - 1 ## making 0 based
    end = int(var_region.split(":")[1].split("-")[1])
    if start >= end:
        raise Exception("start is greater than end in custom position")
    if start < 1:
        raise Exception("invalid start position")
    if end > Chrom().chrom_size(session, chrom):
        raise Exception("invalid end position, greater than chromosome length")  
    if chrom != transcript.chrom:
        raise Exception("custom location is not on the same chromosome as the transcript")
    return (chrom, (start, end))
         
def get_cnv_impact(transcript, var_region, status="new"):
    if var_region in ["CODING", "UTR", "INTRONIC", "GENIC"]:
        regions = transcript.get_requested_region(var_region)
        chrom = transcript.chrom
    else: 
        regions = [parse_region(transcript, var_region)[1]]
        chrom = parse_region(transcript, var_region)[0]
    
    if status == "existing":
        cnv = CNV()
        cnv_list = []
        for region in regions:
            cnv_list = cnv_list + CNV.select_by_location(session, chrom, region[0], region[1])

        if len(cnv_list) == 0:
            print "There is no CNVs in that region please select another variant"
            return False

        print "UID\tSTART\tEND\tVARIANT_TYPE\tCOPY_NUMBER"
        for s in cnv_list:
            print str(s[0].uid) + "\t" + str(s[1].start+1) + "\t" + str(s[1].end) + "\t" + str(s[0].variant_type) + "\t" + str(s[0].copy_number) 

        uid = raw_input("""Input the uid of your CNV: """)
        cnv = CNV.select_by_uid(session, uid)
        if cnv[0].variant_type == "copy_number_loss":
            var_impact = {"CHROM": cnv[1].chrom, "START": cnv[1].start, "END": cnv[1].end, "CNV": -1}
        else:
            var_impact = {"CHROM": cnv[1].chrom, "START": cnv[1].start, "END": cnv[1].end, "CNV": cnv[0].copy_number}
    else:
        print "Available Regions:"
        print "START\tEND"
        for region in regions:
            print str(region[0]+1) + "\t" + str(region[1]) 
        start = int(raw_input("""Input the base 1 start position of your CNV: """))
        end = int(raw_input("""Input the base 1 end position of your CNV: """))
        copy_number = int(raw_input("""Input number of copies you want in your variant where positive numbers are duplications and -1 is a deletion: """))
        if start>=end:
            print "The start position cannot be greater than the end position please start again."
            return False 
        validity = False
        for region in regions: 
            if region[0] <= start and end <= region[1]:
                validity = True
        if validity == False:
            print "The start and end positions must be within a region, please start again."
            return False
        var_impact = {"CHROM": transcript.chrom, "START": start-1, "END": end, "CNV": copy_number}
    return var_impact

def get_str_impact(transcript, var_region, status="new"):
    if var_region in ["CODING", "UTR", "INTRONIC", "GENIC"]:
        regions = transcript.get_requested_region(var_region)
        chrom = transcript.chrom
    else: 
        regions = [parse_region(transcript, var_region)[1]]
        chrom = parse_region(transcript, var_region)[0]
    STR = ShortTandemRepeat()
    str_list = []

    for region in regions:
        str_list = str_list + STR.select_by_location(session, chrom, region[0], region[1])

    new_str_list = []
    if status == "existing":
        for s in str_list: 
            if s[0].pathogenicity != 0: 
                new_str_list.append(s)

    if len(new_str_list) == 0:
        print "There is no STR in that region please select another variant"
        return False

    print "UID\tMOTIF\tPATHOGENICITY"
    for s in new_str_list:
        print str(s[0].uid) + "\t" + str(s[0].motif) + "\t" + str(s[0].pathogenicity) 

    uid = int(raw_input("""Input the uid of your STR: """))
    STR = STR.select_by_uid(session, uid)
    if status != "existing":
        str_length = int(raw_input("""input the size of your short tandem repeat, negative numbers are deletions positive numbers are insersions: """))
    else:
        str_length = STR[0].pathogenicity
    var_impact = {"CHROM": STR[1].chrom, "START": STR[1].start, "END": STR[1].end, "STR": str_length}
    return var_impact

def get_snv_impact(transcript, var_region):
    if var_region == "CODING":
        snv_type = str(raw_input("""Input the impact of your SNV, valid inputs are 'A', 'T', 'G', 'C', 'ANY', 'MISSENSE', 'NONSENSE' or 'SYNONYMOUS': """))
        if snv_type not in ['A', 'T', 'G', 'C', 'ANY', 'MISSENSE', 'NONSENSE', 'SYNONYMOUS']:
            raise ValueError("Not valid impact type for CODING SNV")
    else:
        snv_type = str(raw_input("""Input the impact of your SNV, valid inputs are 'A', 'T', 'G', 'C', or 'ANY': """))
        if snv_type not in ['A', 'T', 'G', 'C', 'ANY']:
            raise ValueError("Not valid impact type for SNV")
    if var_region in ['CODING', 'UTR', 'INTRONIC', 'GENIC']:
        regions = transcript.get_requested_region(var_region)
    else: 
        regions = [parse_region(transcript, var_region)[1]]
    
    location = str(raw_input("""Input the location you want your variant to be.
    The location can be 'ANY' or a location within these ranges 
    inclusive of the start position not inclusive of the end """+ str(regions)+ ": "))   
    if location != "ANY":
        try: 
            location = int(location)
            location_correct = False
            for region in regions: 
                if region[0] <= location < region[1]:
                    location_correct = True
            if location_correct == False:
                raise ValueError("location is not within region")
        except: 
            raise ValueError("location is not a number or 'ANY'")

    return {"SNV_TYPE": snv_type, "LOCATION": location}

def get_indel_impact(transcript, var_region):
    indel_amount = int(raw_input("""input the size of your indel, negative numbers are deletions positive numbers are insersions: """))
    
    if var_region in ['CODING', 'UTR', 'INTRONIC', 'GENIC']:
        regions = transcript.get_requested_region(var_region)
    else: 
        regions = [parse_region(transcript, var_region)[1]]
    
    location = str(raw_input("""Input the location you want your variant to be.
    The location can be 'ANY' or a location within these ranges 
    inclusive of the start position not inclusive of the end """+ str(regions)+ ": "))   
    if location != "ANY":
        try: 
            location = int(location)
            location_correct = False
            for region in regions: 
                if region[0] <= location < region[1]:
                    location_correct = True
            if location_correct == False:
                raise Exception("location is not within region")
        except: 
            raise Exception("location is not a number or 'ANY'")

    return {"INDEL_AMOUNT": indel_amount, "LOCATION": location}

def get_clinvar_impact(transcript, var_region):
    if var_region in ["CODING", "UTR", "INTRONIC", "GENIC"]:
        regions = transcript.get_requested_region(var_region)
        chrom = transcript.chrom
    else: 
        regions = [parse_region(transcript, var_region)[1]]
        chrom = parse_region(transcript, var_region)[0]
    clinvar = ClinVar()
    count = 0 
    clinvar_list = []
    for region in regions:
        clinvar_list = clinvar_list + clinvar.select_by_location(session, chrom, region[0], region[1])
    
    if len(clinvar_list) == 0:
        print "There are no ClinVar variants in that region please select another variant"
        return False
    print "clinvarID\tpos(1-based)\tref\talt\tCLNSIG\tCLNDN"
    for s in clinvar_list:
        clinvar = s[0]
        region = s[1]
        print '\t'.join([str(clinvar.clinvarID), str(region.start+1), str(clinvar.ref), 
        str(clinvar.alt), str(clinvar.CLNSIG), str(clinvar.CLNDN)])

    clinvarID = int(raw_input("""Input the clinvarID of your ClinVar variant: """))
    var_impact = clinvarID
    return {"CLINVAR_ID": var_impact}

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
        var_creation = str(raw_input("Would you like to create a new variant or use an existing variant, input 'new' or 'existing': "))
        if var_creation == "existing":
            var_type = str(raw_input("""What type of variant would you like, options are 'STR', 'ClinVar', 'CNV': """))  # add other variant types
            var_region = str(raw_input("""In what region would you like your variant to be, options are CODING', 'UTR', 'INTRONIC', 'GENIC', or a custom position in the format of chrZ:int-int: """))  # todo add promoter and enhancer
            var_zygosity = str(raw_input("What zygosity would you like your variant to have, options are " + str(zygosity_options) + ": "))
            if var_type == "STR":
                var_impact = get_str_impact(transcript, var_region, "existing")     
            elif var_type == "ClinVar":
                var_impact = get_clinvar_impact(transcript, var_region)
            elif var_type == "CNV": 
                var_impact = get_cnv_impact(transcript, var_region, "existing") 
            else:
                raise ValueError("variant type specified is not valid")
        elif var_creation == "new":
            var_type = str(raw_input("""What type of variant would you like, options are 'SNV', 'INDEL', 'STR', 'CNV': """))  # add other variant types
            var_region = str(raw_input("""In what region would you like your variant to be, options are CODING', 'UTR', 'INTRONIC', 'GENIC', or a custom position in the format of chrZ:int-int: """))  # todo add promoter and enhancer
            var_zygosity = str(raw_input("What zygosity would you like your variant to have, options are " + str(zygosity_options) + ": "))
            if var_type == "INDEL":
                var_impact = get_indel_impact(transcript, var_region)
            elif var_type == "SNV":
                var_impact = get_snv_impact(transcript, var_region)
            elif var_type == "STR":
                var_impact = get_str_impact(transcript, var_region) 
            elif var_type == "CNV": 
                var_impact = get_cnv_impact(transcript, var_region)  
            else:
                raise ValueError("variant type specified is not valid")
        else: 
            raise ValueError("must select new or existing")

    if transcript.chrom in ["chrX", "chrY"] and sex == "XY-MALE" and var_zygosity != "HEMIZYGOUS":
        raise Exception("not valid zygosity.")
    elif var_zygosity not in ["HOMOZYGOUS", "HETEROZYGOUS", "HEMIZYGOUS"]:
        raise Exception("not valid zygosity.")
    var_dict["TYPE"] = var_type
    var_dict["REGION"] = var_region
    var_dict["IMPACT"] = var_impact
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
