# Activate GeneBreaker environment, making the install easier

SIMDIR=$PWD

## Clone the varsim repo, make sure you have git installed
git clone https://github.com/bioinform/varsim

## Run the install
cd $SIMDIR/varsim/
sh build.sh

# Varsim variables
VARSIM_DIR=$PWD/varsim/
OPT_DIR=${VARSIM_DIR}/opt



## Download reference and variant databases 
cd $SIMDIR
### Ref genome version 37
b37_source="http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa"
if [[ ! -f hs37d5.fa ]];then
    wget $b37_source -O  hs37d5.fa
    ${OPT_DIR}/samtools-1.9_install/bin/samtools faidx hs37d5.fa
fi

### Human insertion sequences
if [[ ! -f insert_seq.txt ]];then
    wget http://web.stanford.edu/group/wonglab/varsim/insert_seq.txt
fi

### Supporting large variants
if [[ ! -f GRCh37_hg19_supportingvariants_2013-07-23.txt ]];then
    wget http://web.stanford.edu/group/wonglab/varsim/GRCh37_hg19_supportingvariants_2013-07-23.txt
fi

### dbSNP common variants
if [[ ! -f 00-All.vcf.gz ]];then
    wget -c https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b150_GRCh37p13/VCF/00-common_all.vcf.gz
    gunzip 00-common_all.vcf.gz
fi

############
# OPTIONAL #
############

# Recommended to run with incorporating known variants

# Test varsim install
# Activate environment
source $SIMDIR/varsim/opt/miniconda2/etc/profile.d/conda.sh
conda activate $SIMDIR/varsim/opt/miniconda2/
# sh $SIMDIR/varsim/tests/quickstart_test/quickstart.sh
