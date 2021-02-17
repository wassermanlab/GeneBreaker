# Activate environment
source ../../opt/miniconda3/etc/profile.d/conda.sh
conda activate ../../opt/GeneBreakerEnvironment

# Create subdirectories
mkdir -p Polaris_VCF/
mkdir -p Polaris_VCF/GRCh37
mkdir -p Polaris_VCF/GRCh38

# change to Polaris_VCF directory
cd Polaris_VCF

# Fetch the data
# Download the variants from the online repository
wget -c https://zenodo.org/record/3829960/files/ERR1955420_GRCh37_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR1955420_GRCh37_DeepVariant_0.10.0.vcf.gz.tbi
wget -c https://zenodo.org/record/3829960/files/ERR1955420_GRCh38_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR1955420_GRCh38_DeepVariant_0.10.0.vcf.gz.tbi

wget -c https://zenodo.org/record/3829960/files/ERR1955435_GRCh37_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR1955435_GRCh37_DeepVariant_0.10.0.vcf.gz.tbi
wget -c https://zenodo.org/record/3829960/files/ERR1955435_GRCh38_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR1955435_GRCh38_DeepVariant_0.10.0.vcf.gz.tbi

wget -c https://zenodo.org/record/3829960/files/ERR1955443_GRCh37_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR1955443_GRCh37_DeepVariant_0.10.0.vcf.gz.tbi
wget -c https://zenodo.org/record/3829960/files/ERR1955443_GRCh38_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR1955443_GRCh38_DeepVariant_0.10.0.vcf.gz.tbi

wget -c https://zenodo.org/record/3829960/files/ERR1955495_GRCh37_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR1955495_GRCh37_DeepVariant_0.10.0.vcf.gz.tbi
wget -c https://zenodo.org/record/3829960/files/ERR1955495_GRCh38_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR1955495_GRCh38_DeepVariant_0.10.0.vcf.gz.tbi

wget -c https://zenodo.org/record/3829960/files/ERR1955499_GRCh37_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR1955499_GRCh37_DeepVariant_0.10.0.vcf.gz.tbi
wget -c https://zenodo.org/record/3829960/files/ERR1955499_GRCh38_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR1955499_GRCh38_DeepVariant_0.10.0.vcf.gz.tbi

wget -c https://zenodo.org/record/3829960/files/ERR1955507_GRCh37_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR1955507_GRCh37_DeepVariant_0.10.0.vcf.gz.tbi
wget -c https://zenodo.org/record/3829960/files/ERR1955507_GRCh38_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR1955507_GRCh38_DeepVariant_0.10.0.vcf.gz.tbi

wget -c https://zenodo.org/record/3829960/files/ERR2304556_GRCh37_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR2304556_GRCh37_DeepVariant_0.10.0.vcf.gz.tbi
wget -c https://zenodo.org/record/3829960/files/ERR2304556_GRCh38_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR2304556_GRCh38_DeepVariant_0.10.0.vcf.gz.tbi

wget -c https://zenodo.org/record/3829960/files/ERR2304565_GRCh37_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR2304565_GRCh37_DeepVariant_0.10.0.vcf.gz.tbi
wget -c https://zenodo.org/record/3829960/files/ERR2304565_GRCh38_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR2304565_GRCh38_DeepVariant_0.10.0.vcf.gz.tbi

wget -c https://zenodo.org/record/3829960/files/ERR2304569_GRCh37_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR2304569_GRCh37_DeepVariant_0.10.0.vcf.gz.tbi
wget -c https://zenodo.org/record/3829960/files/ERR2304569_GRCh38_DeepVariant_0.10.0.vcf.gz
wget -c https://zenodo.org/record/3829960/files/ERR2304569_GRCh38_DeepVariant_0.10.0.vcf.gz.tbi


# Remove ref calls
Files=$( ls $PWD/*DeepVariant_0.10.0.vcf.gz )
for file in $Files
do
	echo $file

	# clean it up if it's already there for some reason
	rm ${file::-7}_NoRefCalls.vcf*

	# unzip, and remove RefCall lines 
	gunzip -c $file | grep -v "RefCall" > ${file::-7}_NoRefCalls.vcf 

	# Bgzip and tabix index
	bgzip ${file::-7}_NoRefCalls.vcf 
	tabix -f -p vcf ${file::-7}_NoRefCalls.vcf.gz 

	# move to the right folder
	VCFGZ=${file::-7}_NoRefCalls.vcf.gz
	IFS='/' read -a array <<< $VCFGZ
	IFS='_' read -a array2 <<< ${array[-1]}
	SAMPLE=${array2[0]}
	GENOME=${array2[1]}
	echo $GENOME
	echo $SAMPLE
	mv $VCFGZ $GENOME/${array[-1]}
	mv ${VCFGZ}.tbi $GENOME/${array[-1]}.tbi
done

