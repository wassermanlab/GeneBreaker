import re

columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
##INFO=<ID=MC,Number=.,Type=String,Description="comma separated list of molecular consequence in the form of Sequence Ontology ID|molecular_consequence">
##INFO=<ID=gnomad_exome_af_global,Number=A,Type=Float,Description="Allele Frequency among genotypes, for each ALT allele, in the same order as listed (from /mnt/causes-vnx1/DATABASES/GNOMAD/gnomad.exomes.r2.0.2.sites.norm.vcf.gz)">
##INFO=<ID=gnomad_exome_hom_global,Number=A,Type=Integer,Description="Count of homozygous individuals (from /mnt/causes-vnx1/DATABASES/GNOMAD/gnomad.exomes.r2.0.2.sites.norm.vcf.gz)">
##INFO=<ID=gnomad_genome_af_global,Number=A,Type=Float,Description="Allele Frequency among genotypes, for each ALT allele, in the same order as listed (from /mnt/causes-vnx1/DATABASES/GNOMAD/gnomad.genomes.r2.0.2.sites.wholeGenome.norm.vcf.gz)">
##INFO=<ID=gnomad_genome_hom_global,Number=A,Type=Integer,Description="Count of homozygous individuals (from /mnt/causes-vnx1/DATABASES/GNOMAD/gnomad.genomes.r2.0.2.sites.wholeGenome.norm.vcf.gz)">
##INFO=<ID=CADD,Number=1,Type=String,Description="calculated by self of overlapping values in column 6 from /mnt/causes-vnx1/DATABASES/CADD/whole_genome_SNVs.tsv.gz">

ID = re.compile("\I\D\=[a-xA-Z\_\d\-]*,")
with open("clinvar.norm.anno.vcf") as f:
  for line in f:
    if line.startswith("##INFO="):
      id_str = ID.findall(line)[0][3:-1]
      columns.append(id_str)
new = open("parse.clinvar.anno.tsv", "w+")
new.write("#" + "\t".join(columns)+"\n")

with open("clinvar.norm.anno.vcf") as f:
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
      new.write("\t".join(line_list)+ "\n")
new.close()


# new = open("sample.tsv", "w+")
# new.write("#"+"\t".join(columns))
# with open("clinvar.norm.anno.500.vcf") as f:
#   for line in f:
#     if not line.startswith("#"):
#       line_write = ""
#       fields = line.split("\t")
#       fields[-1] = fields[-1].rstrip()
#       for field in fields[:-1]:
#         line_write = line_write + field + "\t"
        
#       info = field[-1].split(";")
#       info[-1] = field[-1].rstrip()
#       MC="NA"
#       gnomad_exome_af_global="NA"
#       gnomad_exome_hom_global="NA"
#       gnomad_genome_af_global="NA"
#       gnomad_genome_hom_global="NA"
#       CADD = "NA"
#       CLNSIG = "NA"
#       for inf in info:
#         if inf.startswith("MC="):
#           split = inf.split("=")
#           MC = split[1].rstrip()
#         elif inf.startswith("gnomad_exome_af_global="):
#           split = inf.split("=")
#           gnomad_exome_af_global = split[1].rstrip()
#         elif inf.startswith("gnomad_exome_hom_global="):
#           split = inf.split("=")
#           gnomad_exome_hom_global = split[1].rstrip()
#         elif inf.startswith("gnomad_genome_af_global="):
#           split = inf.split("=")
#           gnomad_genome_af_global = split[1].rstrip()
#         elif inf.startswith("gnomad_genome_hom_global="):
#           split = inf.split("=")
#           gnomad_genome_hom_global = split[1].rstrip()
#         elif inf.startswith("CADD="):
#           split = inf.split("=")
#           CADD = split[1].rstrip()
#         elif inf.startswith("CLNSIG="):
#           split = inf.split("=")
#           CLNSIG = split[1].rstrip()
#       line_write = line_write + MC + "\t" + gnomad_exome_af_global + "\t" + gnomad_exome_hom_global + "\t" + gnomad_genome_af_global + "\t" + gnomad_genome_hom_global + "\t" + CADD +"\t" +CLNSIG+"\n"
#       new.write(line_write)
# new.close()
