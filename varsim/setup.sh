#!/bin/bash

module load python/2.7
module load java/1.8.0_121
source activate variant_sim

mkdir varsim_run 
cd varsim_run

set -x

b37_source="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
dbsnp_source="ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz"
varsim_version="0.8.3"
conda=`conda info |grep 'envs directories' |awk '{print $4}'`
bin="${conda}/variant_sim/bin/"

# Download varsim
wget https://github.com/bioinform/varsim/releases/download/v$varsim_version/varsim-$varsim_version.tar.gz
tar xfz varsim-$varsim_version.tar.gz

# Download reference and variant databases 
wget $b37_source -O - | gunzip -c > hs37d5.fa
wget http://web.stanford.edu/group/wonglab/varsim/GRCh37_hg19_supportingvariants_2013-07-23.txt
wget $dbsnp_source -O All.vcf.gz

# index reference
samtools faidx hs37d5.fa