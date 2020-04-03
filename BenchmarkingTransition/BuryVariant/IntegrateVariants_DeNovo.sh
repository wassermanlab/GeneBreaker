# Activate environment
source ../../opt/miniconda3/etc/profile.d/conda.sh
conda activate ../../opt/GeneBreakerEnvironment

# Set up this header 
WORKING_DIR=$PWD
PYTHON=python
REFORMAT_SCRIPT=reformatSimToDeepVariant.py

# Family-based settings

# Proband
IN_PROBAND_GB_VCFA=../../Examples/AD_MSH2_clinvar.vcf
IN_PROBAND_BACKGROUND_VCFGZ=Polaris_VCF/ERR2304597_DeepVariant0.8.noRefCalls.vcf.gz
IN_PROBAND_BACKGROUND_VCF=Polaris_VCF/ERR2304597_DeepVariant0.8.noRefCalls.vcf
PROBAND_COVERAGE=45

# Mother
IN_MOTHER_GB_VCFA=../../Examples/AD_MSH2_clinvar.vcf
IN_MOTHER_BACKGROUND_VCFGZ=Polaris_VCF/ERR1955491_DeepVariant0.8.noRefCalls.vcf.gz
IN_MOTHER_BACKGROUND_VCF=Polaris_VCF/ERR1955491_DeepVariant0.8.noRefCalls.vcf
MOTHER_COVERAGE=45

# Father
IN_FATHER_GB_VCFA=../../Examples/AD_MSH2_clinvar.vcf
IN_FATHER_BACKGROUND_VCFGZ=Polaris_VCF/ERR1955404_DeepVariant0.8.noRefCalls.vcf.gz
IN_FATHER_BACKGROUND_VCF=Polaris_VCF/ERR1955404_DeepVariant0.8.noRefCalls.vcf
FATHER_COVERAGE=45


# Unzip DeepVariant VCFs
gunzip -c $IN_PROBAND_BACKGROUND_VCFGZ > $IN_PROBAND_BACKGROUND_VCF
gunzip -c $IN_FATHER_BACKGROUND_VCFGZ > $IN_FATHER_BACKGROUND_VCF
gunzip -c $IN_MOTHER_BACKGROUND_VCFGZ > $IN_MOTHER_BACKGROUND_VCF


# Modify variant format to match DeepVariant
# This will also output a background file saved as ${IN_PROBAND_BACKGROUND_VCF}_ReformatBackgound.vcf

#PROBAND
$PYTHON $REFORMAT_SCRIPT \
	-I $IN_PROBAND_GB_VCFA \
	-B $IN_PROBAND_BACKGROUND_VCF \
	-C $PROBAND_COVERAGE \
	-F PROBAND \
	-R 

#MOTHER
$PYTHON $REFORMAT_SCRIPT \
        -I $IN_MOTHER_GB_VCFA \
        -B $IN_MOTHER_BACKGROUND_VCF \
        -C $MOTHER_COVERAGE \
        -F MOTHER \
        -R

#FATHER
$PYTHON $REFORMAT_SCRIPT \
        -I $IN_FATHER_GB_VCFA \
        -B $IN_FATHER_BACKGROUND_VCF \
        -C $FATHER_COVERAGE \
        -F FATHER \
        -R

# Organize your files
mv ${IN_PROBAND_BACKGROUND_VCF::-4}_Background_PROBAND.vcf $WORKING_DIR/PROBAND_Background.vcf
mv ${IN_PROBAND_GB_VCFA::-4}_PathoVar_PROBAND.vcf $WORKING_DIR/PROBAND_PathoVar.vcf
mv ${IN_MOTHER_BACKGROUND_VCF::-4}_Background_MOTHER.vcf $WORKING_DIR/MOTHER_Background.vcf
mv ${IN_FATHER_BACKGROUND_VCF::-4}_Background_FATHER.vcf $WORKING_DIR/FATHER_Background.vcf

cd $WORKING_DIR
# BGZIP and tabix the new vcfs
bgzip -f PROBAND_Background.vcf
tabix PROBAND_Background.vcf.gz

bgzip -f PROBAND_PathoVar.vcf
tabix PROBAND_PathoVar.vcf.gz

bgzip -f MOTHER_Background.vcf
tabix MOTHER_Background.vcf.gz

bgzip -f FATHER_Background.vcf
tabix FATHER_Background.vcf.gz

# Concatenate VCFs for relevant individuals
bcftools concat -a PROBAND_Background.vcf.gz PROBAND_PathoVar.vcf.gz -O z -o PROBAND_GeneBreaker.vcf.gz
tabix PROBAND_GeneBreaker.vcf.gz

# Merge across all individuals
bcftools merge -0 PROBAND_GeneBreaker.vcf.gz MOTHER_Background.vcf.gz FATHER_Background.vcf.gz -o Merged_GeneBreaker.vcf

bgzip Merged_GeneBreaker.vcf
tabix Merged_GeneBreaker.vcf.gz

