# Purpose of this script is to take as input a simulated VCF, and a background VCF,
# and then combine the two together.
# Author: Phillip Richmond (phillip.a.richmond@gmail.com)
# Date: 2019-10-18


import sys, argparse


def GetArgs():
	parser = argparse.ArgumentParser()
	Help_detail = "Input Simulated VCF"
        parser.add_argument("-I", "--Input", help=Help_detail, required=True)
	
	Help_detail = "Input Background VCF"
	parser.add_argument("-B", "--Background", help=Help_detail, required=True)

	Help_detail = "Output Combined VCF"
	parser.add_argument("-O", "--Output", help=Help_detail, required=True)

	Help_detail = "Expected Coverage"
	parser.add_argument("-C", "--Coverage", help=Help_detail, required=True, default=40, type=int)

	args = parser.parse_args()
	return args

# Desired output format
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block.">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
# 1       10583   .       G       A       6.1     PASS    .       GT:GQ:DP:AD:VAF:PL      0/1:6:13:10,3:0.230769:4,0,25
def ReadSimVCF(infilename):
	infile = open(infilename,'r')
	for line in infile:
		if line[0]=='#':
			continue
		


def Main():
	ARGS = GetArgs()
	input_VCF_Data = ReadSimVCF(ARGS.Input)
	





