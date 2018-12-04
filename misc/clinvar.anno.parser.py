import re
import argparse
from datetime import date

def parse_args():
    parser = argparse.ArgumentParser(
        description="this script parses clinvar annotated")
    parser.add_argument("-v",  "--vcf_file", help="annotated clinvar")
    args = parser.parse_args()
    return args

def parse_vcf(vcf_file) :
  columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
  ID = re.compile("\I\D\=[a-xA-Z\_\d\-]*,")
  day = str(date.today())
  with open(vcf_file) as f:
    for line in f:
      if line.startswith("##INFO="):
        id_str = ID.findall(line)[0][3:-1]
        columns.append(id_str)
      if line.startswith("##fileDate="):
          day = line.split("=")[1].rstrip()
  annotation_list = ["ANN_Allele", "ANN_Annotation", 
  "ANN_Annotation_Impact", "ANN_Gene_Name", "ANN_Gene_ID", 
  "ANN_Feature_Type", "ANN_Feature_ID"]
  columns = columns + annotation_list

  new = open("parse.clinvar.anno.tsv", "w+")
  new.write("#" + "\t".join(columns)+"\n")

  with open(vcf_file) as f:
    for line in f:
      if not line.startswith("#"):
        new_line = ""
        line_list = ["NA"]*len(columns)
        fields = line.split("\t")
        fields[-1] = fields[-1].rstrip()
        info = fields[-1]
        info = info.split(";")
        info[-1] = info[-1].rstrip()
        ## add fields
        for index,i in enumerate(fields[:-1]):
          line_list[index] = i
        ## add info fields
        for i in info:
          i_split = i.split("=")
          if len(i_split) == 2:
            index = columns.index(i_split[0])
            line_list[index] = i_split[1].rstrip()
            if i_split[0] == "ANN": ## further split the annotation column
              ann = i_split[1].split(",")[0].split("|")
              for idx, val in enumerate(annotation_list):
                index = columns.index(val)
                line_list[index] = ann[idx]
        new.write("\t".join(line_list)+ "\n")
  new.close()

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()
    # Insert ENCODE data to GUD database
    parse_vcf(args.vcf_file)