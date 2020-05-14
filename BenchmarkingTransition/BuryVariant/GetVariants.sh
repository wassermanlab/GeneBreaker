# Activate environment
source ../../opt/miniconda3/etc/profile.d/conda.sh
conda activate ../../opt/GeneBreakerEnvironment

# Fetch the data
mkdir -p Polaris_VCF/
mkdir -p Polaris_VCF/GRCh37
mkdir -p Polaris_VCF/GRCh38

ONLINE_ORIGIN='ZENODO REPOblank'
# Download the variants from the online repository
# wget -c $ONLINE_ORIGIN/sample1.vcf 

# Remove ref calls
Files=$( ls $PWD/Polaris_VCF/*DeepVariant_0.10.0.vcf.gz )
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
	mv $VCFGZ Polaris_VCF/$GENOME/${array[-1]}
	mv ${VCFGZ}.tbi Polaris_VCF/$GENOME/${array[-1]}.tbi
done

