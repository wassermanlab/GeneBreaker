#! /usr/bin/env python

import argparse
import random


def parse_arguments():
    """Parses inputted arguments as described"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-p', '--proband', help="proband's vcf file", type=str, required=True)
    parser.add_argument(
        '-f', '--family',
        help="ped file representing family where the proband is indicated by \
        the sample id = proband", type=str, required=True)
    args = parser.parse_args()
    return args


def output_vcfs(family, vcf):
    "for each family member make a vcf"
    header = ""
    with open(vcf) as f:
        for line in f:
            if line.startswith("##"):
                header = header + line

    for member in family:
        # add header
        # add variants
        member_val = family[member]
        f = open(member + ".vcf", "w+")
        f.write(header)
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + member + "\n")
        if family[member]["var1"] == 0 and family[member]["var2"] == 0:
            None
        elif family[member]["var1"] == family[member]["var2"]:  # homozygous variant
            var = family[member]["var1"][:]
            var[9] = "1/1"
            f.write("\t".join(var))
            f.write("\n")
        elif family[member]["var1"] != 0 and family[member]["var2"] == 0:  # het variant
            var = family[member]["var1"]
            f.write("\t".join(var))
            f.write("\n")
        elif family[member]["var1"] != 0 and family[member]["var2"] != 0:  # het variant
            var = family[member]["var1"]
            f.write("\t".join(var))
            f.write("\n")
            var = family[member]["var2"]
            f.write("\t".join(var))
            f.write("\n")
        f.close()


def make_parents(family):
    "given the proband create the variants for the mom and dad"
    proband = family["Proband"]
    var1 = proband["var1"]
    var2 = proband["var2"]

    if proband['maternal_id'] in family:  # mother is in ped file
        mother = family[proband['maternal_id']]
        if mother["phenotype"] != 1 and mother["phenotype"] != 0:
            mother["var1"] = var1

    if proband['paternal_id'] in family:  # father is in ped file
        father = family[proband['paternal_id']]
        if father["phenotype"] != 1 and father["phenotype"] != 0:
            if var1[0] == "chrY":  # y-linked
                father["var1"] = var1
            elif var2 == 0:  # AD
                father["var1"] = var1
            else:  # AR or AD compound het
                father["var1"] = var2

    return family


def make_siblings(family):
    "given the mom and data create a siblings"
    proband = family["Proband"]
    var1 = proband["var1"]
    var2 = proband["var2"]
    children = family.keys()
    children.remove("Proband")  # remove proband

    if proband['maternal_id'] in family and proband['paternal_id'] in family:  # both parents are in PED
        mother = family[proband['maternal_id']]
        father = family[proband['paternal_id']]
        children.remove(proband['maternal_id'])  # remove mother
        children.remove(proband['paternal_id'])  # remove father
        for child in children:
            if family[child]["phenotype"] == 3:  # child is carrier
                if mother["var1"] != 0 and father["var1"] != 0:  # both parents are carriers
                    family[child]["var1"] = random.choice(
                        (mother["var1"], father["var1"]))
                elif mother["var1"] != 0:  # mom is carrier
                    family[child]["var1"] = mother["var1"]
                elif father["var1"] != 0:  # dad is carrier
                    family[child]["var1"] = father["var1"]

    elif proband['maternal_id'] in family:  # mother is in PED
        mother = family[proband['maternal_id']]
        children.remove(proband['maternal_id'])  # remove mother
        for child in children:
            if family[child]["phenotype"] == 3:  # child is carrier
                family[child]["var1"] = mother["var1"]

    elif proband['paternal_id'] in family:  # father is in PED
        father = family[proband['paternal_id']]
        children.remove(proband['paternal_id'])  # remove father
        for child in children:
            if family[child]["phenotype"] == 3:  # child is carrier
                family[child]["var1"] = father["var1"]

    else:  # no parent is in the family
        for child in children:
            if family[child]["phenotype"] == 3:  # child is carrier
                family[child]["var1"] = random.choice((var1, var2))

    for child in children:
        if family[child]["phenotype"] == 2 and var1[0] != "chrX":  # child is affected
            family[child]["var1"] = var1
            family[child]["var2"] = var2
        if family[child]["phenotype"] == 2 and var1[0] == "chrX":
            if family[child]["sex"] == 2:  # child is female
                family[child]["var1"] = var1
                family[child]["var2"] = var2
            if family[child]["sex"] == 1:  # child is male
                # mom always gets var1 anyways so its matching :D
                family[child]["var1"] = var1

    return family


def main():
    args = parse_arguments()
    proband_vcf = args.proband
    family_ped = args.family
    var1 = 0
    var2 = 0
    proband_vcf = [line.rstrip() for line in open(proband_vcf)]
    for line in proband_vcf:
        if not line.startswith("#"):
            ln = line.split("\t")
            if var1 == 0:
                if ln[9] == "0/1":
                    var1 = ln
                elif ln[9] == "1/1":
                    ln[9] = "0/1"
                    var1 = ln
                    var2 = ln
                elif ln[9] == "1":
                    var1 = ln
            elif var2 == 0:
                var2 = ln
            else:
                raise Exception("too many variants")
    family = {}
    family_ped = [line.rstrip() for line in open(family_ped)]
    for line in family_ped:
        if not line.startswith("#"):
            ln = line.split("\t")
            family[ln[1]] = {"paternal_id": ln[2],
                             "maternal_id": ln[3],
                             "sex": int(ln[4]),
                             "phenotype": int(ln[5]),
                             "var1": 0,
                             "var2": 0}
    if "Proband" in family:
        family["Proband"]["var1"] = var1
        family["Proband"]["var2"] = var2
        family = make_parents(family)
        family = make_siblings(family)
        output_vcfs(family, args.proband)
    else:
        raise Exception("no 'Proband' sample_id in PED file")

if __name__ == "__main__":
    main()
