# Purpose of this script is to take as input a simulated VCF, and produces a bamsurgeon varfile.
# Author: Phillip Richmond (phillip.a.richmond@gmail.com)


import sys, argparse


def GetArgs():
    parser = argparse.ArgumentParser()
    Help_detail_Input = "Input Simulated VCF(s)"
    Help_detail_Output = "Output Bamsurgeon varfile"
    parser.add_argument("-I", "--Input", help=Help_detail_Input, required=True)
    parser.add_argument("-O", "--Output", help=Help_detail_Output, required=True)
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
def ReadSimVCFWriteBamSurgeon(infilename,outfilename):
    infile = open(infilename,'r')
    outfile = open(outfilename,'w')
    for line in infile:
        if line[0]=='#':
            continue
        cols = line.strip('\n').split('\t')
        chrom = cols[0]
        # a catch for empty lines
        try:
            pos = int(cols[1])
        except(IndexError):
            break
        identifier = cols[2]
        ref = cols[3]
        alt = cols[4]
        genotype = cols[9].split(':')[0]
        if genotype=='0/1':
            ratio=0.5
        elif genotype=='1/1':
            ratio=1.0
        elif (genotype=='0/0') or (genotype=='./.'):
            print("No variant needed for simulation")
            continue
       # if indel
        # output looks like
        # chrom start   end ratio   INS||DEL    insertedSequence
        # start and end are 0-based
        # so start=pos-1, and end=(pos-1)+len(ref) for deletions
        # and
        # start=pos-1, end=pos for insertions
        start=pos-1
        if len(ref) > len(alt):
            Type='DEL'
            end=start+(len(ref)-len(alt))
            outfile.write("%s\t%d\t%d\t%f\t%s\t%s\n"%(chrom,start,end,ratio,Type,alt))
 
        elif len(ref) < len(alt):
            Type='INS'
            end=pos
            outfile.write("%s\t%d\t%d\t%f\t%s\t%s\n"%(chrom,start,end,ratio,Type,alt))
 
        elif len(ref) == len(alt):
            # if snv
        # output looks like:
        # pos is 1-based here
        # chrom pos pos ratio    alt_allele
            outfile.write("%s\t%d\t%d\t%.1f\t%s\n"%(chrom,pos,pos,ratio,alt))




def Main():
    ARGS = GetArgs()
    ReadSimVCFWriteBamSurgeon(ARGS.Input,ARGS.Output)

if __name__=="__main__":
    Main()
