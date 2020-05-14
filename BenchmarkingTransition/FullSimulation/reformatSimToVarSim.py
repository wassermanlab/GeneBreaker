
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
