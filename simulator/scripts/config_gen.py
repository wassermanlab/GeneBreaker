# todo 
# fix TRIO responses
import argparse
import json
from GUD.ORM import Gene
import os, sys 
from sqlalchemy import create_engine, Index
from sqlalchemy.orm import Session
import warnings

# Establish a SQLalchemy session w/ GUD
db_name = "mysql://{}:@{}:{}/{}".format("ontarget_r",
    "ontarget.cmmt.ubc.ca", "5506", "hg19")
try:
    engine = create_engine(db_name, echo=False)
    session = Session(engine)
except:
    raise ValueError("Cannot connect to GUD: %s" % "hg19")

def get_variant_dict():
    var_dict = dict()
    var_type = str(raw_input("""What type of variant would you like, options are
    'SNV', 'INDEL': """)) ## add other variant types
    var_region = str(raw_input("""In what region would you like your variant to be, options are
    CODING', 'UTR', 'INTRONIC': """))  ## todo add promoter and enhancer
    if var_type == "INDEL":
        var_impact = int(raw_input("""input the size of your indel, negative
        numbers are deletions positive numbers are insersions: """))  
    elif var_type == "SNV":
        if var_region == "CODING":
            var_impact = str(raw_input("""Input the impact of your SNV, valid
            inputs are 'A', 'T', 'G', 'C', 'ANY', 'MISSENSE', 'NONSENSE' or 'SILENT': """))
            if var_impact not in ['A', 'T', 'G', 'C', 'ANY', 'MISSENSE', 'NONSENSE', 'SILENT']:
                raise Exception("Not valid impact type for CODING SNV")
        else: 
            var_impact = str(raw_input("""Input the impact of your SNV, valid
            inputs are 'A', 'T', 'G', 'C', or 'ANY': """))
            if var_impact not in ['A', 'T', 'G', 'C', 'ANY']:
                raise Exception("Not valid impact type for SNV")
    else:
        raise Exception("variant type specified is not valid")
    var_location = str(raw_input("Specify a location or input 'ANY': "))
    if var_location == 'ANY':
        None 
    elif var_location.isdigit():
        var_location = int(var_location)
    else: 
        raise Exception("Not valid variant location")
    var_dict["TYPE"] = var_type
    var_dict["REGION"] = var_region
    var_dict["IMPACT"] = {"TYPE_IMPACT": var_impact, "LOCATION": var_location}
    return var_dict


def get_main_dict():
    config_dict = dict()
    ## highly coupled with GUD, need to figure out a better way to do this 
    gene_sym = str(raw_input("Enter the gene symbol of the gene you would like to use: "))
    gene = Gene()
    genes = gene.select_by_name(session, gene_sym) ## GALK2
    if len(genes) == 0:
        raise Exception("No genes by that gene symbol")
    counter = 0
    for g in genes: 
        print str(counter) + "\t" + g.name + "\t" + str(g.txStart) + "\t" + str(g.txEnd)
        counter = counter + 1
    gene_index = int(raw_input("what is the index of the transcript would you like to use: "))
    gene = genes[gene_index] ## get gene that we need 
    config_dict["GENE_NAME"] = gene.name2
    config_dict["CHR"] = gene.chrom
    config_dict["STRAND"] = gene.strand
    config_dict["TXSTART"] = gene.txStart
    config_dict["TXEND"] = gene.txEnd
    inheritance = str(raw_input("""What inheritance model would you like to use?            
    valid types are 'DE-NOVO', 'BI-PARENTAL', 'MATERNAL', 'PATERNAL': """))
    trio = str(raw_input("""Would you like to output a trio or just the child, 
    valid types are 'TRIO' or 'SINGLE': """))
    if inheritance in ['DE-NOVO', 'BI-PARENTAL', 'MATERNAL', 'PATERNAL']:
        config_dict["INHERITANCE"] = inheritance
    else:
        raise Exception("Inheritance is not in ['DE-NOVO', 'BI-PARENTAL', 'MATERNAL', 'PATERNAL']")
    if trio in ['TRIO', 'SINGLE']:
        config_dict["TRIO"] = trio
    else:
        raise Exception("Trio entered is not in ['TRIO', 'SINGLE']")
    return config_dict

def main():
    config_name = str(raw_input("enter the name of your config file: "))
    config_dict = get_main_dict()
    print("""===========Variant I===========""")
    # get variant I 
    config_dict["VAR1"] = get_variant_dict()
    print("""===========Variant II===========""")
    # get variant II
    var2_ans = str(raw_input("Would you like a second variant [y/n]: "))
    if var2_ans not in  ["y", "n"]:
        raise Exception("Not valid response.") 
    elif var2_ans == "y":
        config_dict["VAR2"] = get_variant_dict()
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