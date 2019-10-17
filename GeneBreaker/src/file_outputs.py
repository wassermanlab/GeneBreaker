import json
import datetime

# both take in a json file like this {var1: {}, var2: {}, family: {}}
# vars have the following keys: chrom, pos, id, ref, alt, qual, filter, info, format, proband
# family is like this {member_id: {relationship, sex, var1, var2, affected}}
# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype (unaffected=1; affected=2)

def get_row(individual_id, individual, mother, father):
    row = "FAM\t" + individual_id + "\t"
    if (individual_id in ["father", "mother"]):
        row = row + "0\t0\t"
    else:
        if (father):
            row = row + "father\t"
        else:
            row = row + "0\t"
        if (mother):
            row = row + "mother\t"
        else:
            row = row + "0\t"
    if (individual["sex"] == "XX"):
        row = row + "2\t"
    else:
        row = row + "1\t"
    if (individual["affected"] == True):
        row = row + "2\n"
    else:
        row = row + "1\n"
    return row


def make_ped(variants, filename):
    header = "#Family ID\tIndividual ID\tPaternal ID\tMaternal ID\tSex\tPhenotype\n"
    father = 'father' in variants['family']
    mother = 'mother' in variants['family']
    rows = ""
    for key, val in variants['family'].items():
        rows = rows + get_row(key, val, mother, father)
    f = open(filename, "w+")
    f.write(header+rows)
    f.close()


def make_vcf(variants, filename):
    header = "##fileformat=VCFv4.2\n"
    header = header + "##fileDate=" + str(datetime.date.today()) + "\n"
    header = header + "##source=variant_simulator\n"
    header = header + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n"
    header = header + "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n"
    header = header + "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
    header = header + "##ALT=<ID=DEL,Description=\"Deletion\">\n"
    header = header + "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">\n"
    header = header + "##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">\n"
    header = header + "##ALT=<ID=INS:ME:LINE,Description=\"Insertion of LINE1 element\">\n"
    header = header + "##ALT=<ID=INS:ME:SVA,Description=\"Insertion of SVA element\">\n"
    header = header + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\"\n"
    header = header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    row1 = "\t".join([variants["var1"]["chrom"], variants["var1"]["pos"], variants["var1"]["id"],
                      variants["var1"]["ref"], variants["var1"]["alt"], variants["var1"]["qual"],
                      variants["var1"]["filter"], variants["var1"]["info"], variants["var1"]["format"]])
    row2 = ""
    chrom = variants['var1']['chrom']
    if (variants["var2"] != ""):  # no second variant ie. 1 variant or 1 homozygous/hemizygous variant
        row2 = "\t".join([variants["var2"]["chrom"], variants["var2"]["pos"], variants["var2"]["id"],
                          variants["var2"]["ref"], variants["var2"]["alt"], variants["var2"]["qual"],
                          variants["var2"]["filter"], variants["var2"]["info"], variants["var2"]["format"]])
    
    for key, val in variants['family'].items():
        header = header + "\t" + key
        if val["sex"] == "XY" and chrom in ["chrX", "chrY"]:
            if val["var1"]:
                row1 = row1 + "\t1/1"
            elif not val["var1"]:
                row1 = row1 + "\t0/0"
            if val["var2"] and row2 == "":
                row1 = row1 + "\t1/1"
            elif row2 == "":
                row1 = row1 + "\t0/0"
            if val["var2"] and row2 != "":
                row2 = row2 + "\t1/1"
            elif row2 != "":
                row2 = row2 + "\t0/0"
        else: 
            if val["var1"] and val["var2"] and row2 == "": ## homozygous on var 1 
                row1 = row1 + "\t1/1"
            elif (val["var1"] or val["var2"]) and row2 == "": ## heterozygous on var 1 
                row1 = row1 + "\t0/1"
            elif not val["var1"] and not val["var2"] and row2 == "": ## homozygous on var 1 
                row1 = row1 + "\t0/0"

            if val["var1"] and row2 != "": ## het on var 1 
                row1 = row1 + "\t0/1"
            elif not val["var1"] and row2 != "": ## none var 1
                row1 = row1 + "\t0/0"
            if val["var2"] and row2 != "": ## het on var 2
                row2 = row2 + "\t0/1"
            elif not val["var2"] and row2 != "": ## none var 2
                row2 = row2 + "\t0/0"

    f = open(filename, "w+")
    f.write(header + "\n" + row1 + "\n"+row2 + "\n")
    f.close()

    #     header = header + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\"\n"
    #     header = header + "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\"\n"
    #     header = header + "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\"\n"
    #     header = header + "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\"\n"
    #     header = header + "##ALT=<ID=DUP,Description=\"Duplication\"\n"
    #     header = header + "##ALT=<ID=DEL,Description=\"Deletion\"\n"


variants = {'var1': {'alt': 'T',
                     'chrom': 'chr17',
                     'filter': '.',
                     'format': 'GT',
                     'id': '450265',
                     'info': '.',
                     'pos': '72121409',
                     'proband': '1/1',
                     'qual': '.',
                     'ref': 'C'},
            'var2': '',
            'family': {'mother': {'relationship': 'mother', 'sex': 'XX', 'var1': True, 'var2': False, 'affected': True},
                       'father': {'relationship': 'father', 'sex': 'XY', 'var1': False, 'var2': True, 'affected': False},
                       'sister1': {'relationship': 'sibling', 'sex': 'XX', 'var1': False, 'var2': False, 'affected': False},
                       'proband': {'relationship': 'sibling', 'sex': 'XX', 'var1': True, 'var2': True, 'affected': True}}}
variants = {'var1': {'alt': 'T',
                     'chrom': 'chrX',
                     'filter': '.',
                     'format': 'GT',
                     'id': '450265',
                     'info': '.',
                     'pos': '72121409',
                     'proband': '1/1',
                     'qual': '.',
                     'ref': 'C'},
            'var2': '',
            'family': {'mother': {'relationship': 'mother', 'sex': 'XX', 'var1': True, 'var2': False, 'affected': True},
                       'father': {'relationship': 'father', 'sex': 'XY', 'var1': False, 'var2': True, 'affected': False},
                       'sister1': {'relationship': 'sibling', 'sex': 'XX', 'var1': False, 'var2': False, 'affected': False},
                       'proband': {'relationship': 'sibling', 'sex': 'XX', 'var1': True, 'var2': True, 'affected': True}}}
variants = {'var1': {'alt': 'T',
                     'chrom': 'chrX',
                     'filter': '.',
                     'format': 'GT',
                     'id': '450265',
                     'info': '.',
                     'pos': '72121409',
                     'proband': '0/1',
                     'qual': '.',
                     'ref': 'C'},
            'var2': {'alt': 'A',
                     'chrom': 'chrX',
                     'filter': '.',
                     'format': 'GT',
                     'id': '450265',
                     'info': '.',
                     'pos': '72121409',
                     'proband': '0/1',
                     'qual': '.',
                     'ref': 'C'},
            'family': {'mother': {'relationship': 'mother', 'sex': 'XX', 'var1': True, 'var2': False, 'affected': True},
                       'father': {'relationship': 'father', 'sex': 'XY', 'var1': False, 'var2': True, 'affected': False},
                       'sister1': {'relationship': 'sibling', 'sex': 'XX', 'var1': False, 'var2': False, 'affected': False},
                       'proband': {'relationship': 'sibling', 'sex': 'XX', 'var1': True, 'var2': True, 'affected': True}}}

make_ped(variants, 'test.vcf')
