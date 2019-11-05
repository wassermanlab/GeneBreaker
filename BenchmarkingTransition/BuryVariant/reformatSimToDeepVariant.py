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



def ReadSimVCF(infilename):
	infile = open(infilename,'r')
	for line in infile:
		if line[0]=='#':
			continue
		chrom,


def Main():
	ARGS = GetArgs()
	input_VCF_Data = ReadSimVCF(ARGS.Input)
	





