# Purpose of this script is to take as input a simulated VCF, and a background VCF,
# and then combine the two together.
# Author: Phillip Richmond (phillip.a.richmond@gmail.com)
# Date: 2019-10-18


import sys, argparse
from math import floor
from math import ceil



def GetArgs():
	parser = argparse.ArgumentParser()
	Help_detail = "Input Simulated VCF(s)"
	parser.add_argument("-I", "--Input", help=Help_detail, required=True)
	
	Help_detail = "Input Background VCF"
	parser.add_argument("-B", "--Background", help=Help_detail, required=True)

	Help_detail = "Expected Coverage"
	parser.add_argument("-C", "--Coverage", help=Help_detail, required=True, default=45, type=int)

	Help_detail = "Chr Removal, set this option if you want to remove the 'chr' from chromosome names"
	parser.add_argument("-R", "--RemoveChr", help=Help_detail, action="store_true")

	Help_detail = "Family Member, choose one of [PROBAND, MOTHER, FATHER]"
	parser.add_argument("-F", "--FamilyMember", help=Help_detail, type=str, required=True)
	args = parser.parse_args()
	return args


# Current output from GeneBreaker
##fileformat=VCFv4.2
##fileDate=20090805
##source=variant_simulator
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	PROBAND
# chr2	47641508	snv_47641508_MISSENSE	A	C	.	.	.	GT	0/1 

# Here We will read in the simulator VCF and output info on the varant(s)
def ReadSimVCF(infilename):
	infile = open(infilename,'r')
	Vars=[]
	for line in infile:
		if line[0]=='#':
			continue
		cols=line.strip('\n').split('\t')
		chrom=cols[0]
		pos = cols[1]
		identifier = cols[2]
		ref = cols[3]
		alt = cols[4]
		genotype = cols[9]
		var_id = "%s;%s;%s;%s;%s;%s"%(chrom,pos,identifier,ref,alt,genotype)
		Vars.append(var_id)
	print("Read in these variants:")
	print("chrom;pos;identifier;ref;alt;genotype")
	for each in Vars:
		print(each)
	return Vars


# This function will take in a variant with positional and genotype information,
# and output a formatted VCF line. Examples below

# Desired output format
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block.">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
# 1	   10583   .	   G	   A	   6.1	 PASS	.	   GT:GQ:DP:AD:VAF:PL	  0/1:6:13:10,3:0.230769:4,0,25
# Some good quality ones
# 1	13824118	.	G	A	59.7	PASS	.	GT:GQ:DP:AD:VAF:PL	0/1:60:44:18,26:0.590909:59,0,79
# 1	13829419	.	G	A	50.3	PASS	.	GT:GQ:DP:AD:VAF:PL	0/1:50:43:21,22:0.511628:50,0,150
# 1	13824080	.	C	CAA	38.7	PASS	.	GT:GQ:DP:AD:VAF:PL	0/1:39:41:18,21:0.512195:38,0,55
# Good quality homo
# 1	13844997	.	C	T	150	PASS	.	GT:GQ:DP:AD:VAF:PL	1/1:66:32:0,32:1:150,65,0

def ReformatVarLine(variant_id,remove_chr,coverage):
	sim_chrom,sim_pos,sim_identifier,sim_ref,sim_alt,sim_genotype = variant_id.split(';')
	if (remove_chr):
		new_sim_chrom = sim_chrom[3:]
	else:
		new_sim_chrom = sim_chrom
	# These are fields for the new output
	CHROM = new_sim_chrom
	POS = sim_pos
	ID = sim_identifier
	REF = sim_ref
	ALT = sim_alt
	QUAL = '100'
	GT = sim_genotype
	GQ = '75'
	DP = str(coverage)
	# for het genotypes, with synthetic PL scores, make the alt depth higher if odd coverage using floor/ceil
	if GT == '0/1':
		ref_depth = floor(coverage/2)
		alt_depth = ceil(coverage/2)
		PL = '100,0,100'
	# for homo genotypes, ref depth 0
	if GT == '1/1':
		ref_depth = 0
		alt_depth = coverage
		PL = '150,100,0'
	AD = '%d,%d'%(ref_depth,alt_depth)
	VAF = float(alt_depth/coverage)
	newline = "%s\t%s\t%s\t%s\t%s\t%s\tPASS\t.\tGT:GQ:DP:AD:VAF:PL\t%s:%s:%s:%s:%0.3f:%s\n"%(CHROM, POS, ID, REF, ALT, QUAL, GT, GQ, DP, AD, VAF, PL) 
	print(newline)
	return(newline)
	
	


# This function will reformat a background VCF, and the Simulated VCFs to match
# Then the reformatted VCF will be ready for concatenation
def ReformatVCFs(simulated_variants,background_vcf, infilename, family_member, remove_chr, coverage):
	# The output files will be the same as input but add family member and PathoVar or Background 
	backgroundvcf = open(background_vcf,'r')
	out_backgroundvcf_filename = "%s_Background_%s.vcf"%(background_vcf[:-4],family_member)
	out_backgroundvcf = open(out_backgroundvcf_filename,'w')
	outfilename = "%s_PathoVar_%s.vcf"%(infilename[:-4],family_member)
	outfile = open(outfilename,'w')
	print("Writing to these outfiles: ")
	print(outfilename)
	print(out_backgroundvcf_filename)

		
	# Parse through header and print to outfiles
	# Then parse through rest of file and make a full background VCF labeled with family member
	for line in backgroundvcf:
		# ignore header info, just write to outfiles
		if line[0:2]=='##':
			out_backgroundvcf.write(line)
			outfile.write(line)
			continue
		# Re-make the label here to match the family member
		if line[0]=='#':
			cols=line.strip('\n').split('\t')
			cols[-1] = family_member 
			newline="\t".join(cols)
			out_backgroundvcf.write("%s\n"%newline)
			outfile.write("%s\n"%newline)
			continue
		out_backgroundvcf.write(line)	

	# Reformat the line with the function above
	# Then print to outfile
	for var in simulated_variants:
		newline = ReformatVarLine(var,remove_chr,coverage)
		outfile.write(newline)
		

def Main():
	ARGS = GetArgs()
	# Account for multiple VCFs
	# Create a list where I store variants from each VCF (allowing each VCF multiple vars)
	input_VCF_Variants = ReadSimVCF(ARGS.Input)
	ReformatVCFs(input_VCF_Variants, ARGS.Background, ARGS.Input, ARGS.FamilyMember, ARGS.RemoveChr, ARGS.Coverage)

if __name__=="__main__":
	Main()
