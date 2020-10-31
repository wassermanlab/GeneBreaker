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

	Help_detail = "Family Member, choose one of [proband, mother, father]"
	parser.add_argument("-F", "--FamilyMember", help=Help_detail, type=str, required=True)
	args = parser.parse_args()
	return args


# Current output from GeneBreaker
# 
# ##fileformat=VCFv4.2
# ##fileDate=20090805
# ##source=variant_simulator
# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	proband
# X	31200303	CNV_31200302	C	<DUP:TANDEM>	.	.	SVTYPE=DUP;END=31226059;SVLEN=25758;	GT	0/1
# 

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
                svinfo = cols[7]
		genotype = cols[9]
		var_id = "%s;%s;%s;%s;%s;%s;%s"%(chrom,pos,identifier,ref,alt,svinfo,genotype)
		Vars.append(var_id)
	print("Read in these variants:")
	print("chrom;pos;identifier;ref;alt;svinfo;genotype")
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
# - This is a listed insertion seq - #
# chr19	3145577	polaris_manta:chr19_3145577_0_498_MantaINS:79241:0:0:0:0:0	T	TCTCACTCTGTCACCCAGGCTGGAGTGCAGTGGTGTGATCTCGGCTCACTGCAACCTCTGCCTCCCAGGTTGAAGCAATTTTCCTGCCTCAGCCTCCCGAGTAGCTGGGATTACTGGCGCCCACGACCACGCCCAGCTAATTTTTGTATTTTTAGTAGTAATTTTTGTATTTCACCATATTGGTCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATCCACCCGCCTTGGCCTCCCAAAGTGCTGGGATTACAGGCATGAGCCACTGTGCCTGGCCCCTGATACCTTTTCCTCCCTAAATACTCAGCGAGTTTGAGTTGCACTGAGAAGGAAAAACACTGCGGGACCTTGATCAAAAGCATGACACAGTGCCCTCCGCTAACACACAGTCCTGATTTGCAATTTGCCAGTTTTCCCAATTTCCGTGTGAATGTCCTTTACAGCAAGAGAAAAGCCTTTTTCTTTATTCCTTTTTTTTTTTTTTTTTTTGGAGACGGAGC	.	PASS	SVTYPE=INS;SVLEN=498;SOURCE=polaris_manta;MSOURCE=polaris_manta;PEDIGREE=inconsistent;HAMMING=1;TRIO_CALL_RATE=1;NUM_ALT_KIDS=50;FRAC_INCONSI=0.02;AC=2;AN=2	GT	1/1
#
# - This is a listed DEL - # 
# chr19	3173223	MantaDEL:206077:1:2:0:0:0_78	T	<DEL>	.	PEDIGREE;PASS_RATE;HWE;INFERIOR	END=3175784;SVTYPE=DEL;SVLEN=-2561;SOURCE=mantaPG;MSOURCE=mantaPG;PEDIGREE=inconsistent;HAMMING=2;TRIO_CALL_RATE=1;NUM_ALT_KIDS=25;FRAC_INCONSI=0.04;CALL_RATE=0.92;PASS_RATE=0.35;HWE=0;ALT_AF=0.323;AC=1;AN=2	GT	0/1
#
#

# Purpose of this script is to take as input a simulated VCF, and produces a VarSim-formatted VCF.
# Only real change right now is for MEIs, as I want the actual insert sequence
# Author: Phillip Richmond (phillip.a.richmond@gmail.com)
# Date: 2020-05-04
# May the fourth be with you
# I should be writing my thesis right now


import sys, argparse
from Bio import SeqIO


def GetArgs():
	parser = argparse.ArgumentParser()
	Help_detail = "Input Simulated VCF(s)"
	parser.add_argument("-I", "--Input", help=Help_detail, required=True)
	
	args = parser.parse_args()
	return args


# Current output from GeneBreaker
##fileformat=VCFv4.2
##fileDate=20090805
##source=variant_simulator
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	proband
# chr2	47641508	snv_47641508_MISSENSE	A	C	.	.	.	GT	0/1 

# Here We will read in the simulator VCF and output info on the varant(s)
def ReadSimVCF(infilename,mei_fasta_file):
	infile = open(infilename,'r')
	outfile = open("%s_forVarSim.vcf"%(infilename[:-4]),'w')
	for line in infile:
		if line[0]=='#':
			outfile.write(line)
			continue
		cols = line.strip('\n').split('\t')
		chrom = cols[0]
		# a catch for empty lines
		try:
			pos = cols[1]
		except(IndexError):
			break
		identifier = cols[2]
		ref = cols[3]
		alt = cols[4]
		genotype = cols[9]

		# Catch for MEIs
		if 'mei' in identifier:
			varsim_mei_line = ReformatMEILine(line,mei_fasta_file)
			outfile.write(varsim_mei_line)
			continue
		else:
			outfile.write(line)


# Original line looks like this:
# 6	74320297	mei_74320297_LINE	A	<INS:MEI:LINE>	.	.	SVTYPE=INS;END=74326316;SVLEN=6019;	GT	1/1

# New line needs to look like this:
# 6	74320297	mei_74320297_LINE	A	ACCGTC...LINESequence...	.	.	.	GT	0/1
def ReformatMEILine(line,mei_fasta_file):
	# collect seqs into dictionary from mei fasta file
	record_dict = SeqIO.to_dict(SeqIO.parse(mei_fasta_file, "fasta"))

	# parse the line
	cols = line.strip('\n').split('\t')
	# no changes here
	chrom=cols[0]
	pos = cols[1]
	identifier = cols[2]
	ref = cols[3]
	genotype = cols[9]
	
	# now for a new alt
	mei_type = identifier.split('_')[-1]
	new_alt = record_dict[mei_type].seq
	
	# print out	
	newline = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(chrom, pos, identifier, ref, new_alt, '.', '.', '.', 'GT', genotype)
	return newline

def Main():
	ARGS = GetArgs()
	MEI_FASTA='./MEI_Seqs/MEIs.fa'
	ReadSimVCF(ARGS.Input,MEI_FASTA)

if __name__=="__main__":
	Main()









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
	outfilename = "%s.vcf"%(infilename[:-4])
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
