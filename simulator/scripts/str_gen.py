#! /usr/bin/env python

import argparse

def parse_arguments():
    """Parses inputted arguments as described"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s', '--str_bed', help="bedfile of strs", type=str)
    parser.add_argument(
        '-o', '--output', help='type of output together yields one vcf seperate yields multiple', 
        type=str, default="together", choices=['together', 'seperate'])
    parser.add_argument(
        '-n', '--name', help='name of output vcf', type=str, default="str_output")
    args = parser.parse_args()
    return args

def get_retraction(line):
    """input is line array [chrom, start, stop, len, motif, full motif, exp] 
    from bed file. ouput is array [chrom, pos, ID, ref, alt]
    """
    chrom = line[0]
    pos = str(int(line[1]) + 1)
    ID = chrom + "_" + pos + "_" + line[6]
    multiple = -1*int(line[6])
    if multiple*int*line[3]>len(line[5]):
        raise Exception("requested retraction is larger than short tandem repeat" )
    ref = line[4]*multiple
    alt = ""
    return [chrom, pos, ID, ref, alt]

def get_expantion(line):
    """input is line array [chrom, start, stop, len, motif, full motif, exp] 
    from bed file. ouput is array [chrom, pos, ID, ref, alt]
    """
    chrom = line[0]
    pos = str(int(line[1]) + 1)
    ID = chrom + "_" + pos + "_" + line[6]
    ref = line[4][0]
    alt = line[4]*int(line[6])+ line[4][0]
    return [chrom, pos, ID, ref, alt]

def parse_rows(str_bed):
    """input is the bed file with the following columns:
    chrom, start, stop, len, motif, full motif, exp.
    ouput is array or line arrays for the vcf."""
    vcf = []
    with open(str_bed) as f:
        for line in f:
            if not line.startswith("#"):
                split_line = line.split("\t")
                split_line[len(split_line)-1] = split_line[len(split_line)-1].rstrip()
                if split_line[len(split_line)-1] < 0:
                    vcf.append(get_retraction(split_line))
                elif split_line[len(split_line)-1] > 0:
                    vcf.append(get_expantion(split_line))
    return vcf

def output_rows(output, vcf_rows, name):
    """input is the type of output. output is one or many vcfs."""
    if output is "together":
        h = "#CHROM\tPOS\tID\tREF\tALT\n"
        f = open(name+".vcf","w+")
        f.write(h)
        for row in vcf_rows: 
            f.write("\t".join(row))
            f.write("\n")
        f.close()
    else: 
        counter = 1 
        for row in vcf_rows: 
            f = open(str(counter)+name+".vcf","w+")
            f.write("#CHROM\tPOS\tID\tREF\tALT\n")
            f.write("\t".join(row))
            f.write("\n")
            f.close()
            counter = counter + 1

def main():
    args = parse_arguments()
    str_bed = args.str_bed
    output = args.output
    name = args.name
    vcf_rows = parse_rows(str_bed)
    output_rows(output, vcf_rows, name)

if __name__ == "__main__":
    main()