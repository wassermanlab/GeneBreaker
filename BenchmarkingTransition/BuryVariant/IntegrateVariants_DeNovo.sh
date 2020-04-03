# Activate environment
source ../../opt/miniconda3/etc/profile.d/conda.sh
conda activate ../../opt/GeneBreakerEnvironment

# Set up this header 
WORKING_DIR=$PWD
PYTHON=python
REFORMAT_SCRIPT=reformatSimToDeepVariant.py

# Family-based settings
# Merged VCF file
MERGED_GB_VCF=../../Examples/MSH2_AutosomalDominant_DeNovo.vcf

# Proband
IN_PROBAND_GB_VCFA=${MERGED_GB_VCF::-4}.proband.vcf
IN_PROBAND_BACKGROUND_VCFGZ=Polaris_VCF/ERR2304597_DeepVariant0.8.noRefCalls.vcf.gz
IN_PROBAND_BACKGROUND_VCF=Polaris_VCF/ERR2304597_DeepVariant0.8.noRefCalls.vcf
PROBAND_COVERAGE=45

# Mother
IN_MOTHER_GB_VCFA=${MERGED_GB_VCF::-4}.mother.vcf
IN_MOTHER_BACKGROUND_VCFGZ=Polaris_VCF/ERR1955491_DeepVariant0.8.noRefCalls.vcf.gz
IN_MOTHER_BACKGROUND_VCF=Polaris_VCF/ERR1955491_DeepVariant0.8.noRefCalls.vcf
MOTHER_COVERAGE=45

# Father
IN_FATHER_GB_VCFA=${MERGED_GB_VCF::-4}.father.vcf
IN_FATHER_BACKGROUND_VCFGZ=Polaris_VCF/ERR1955404_DeepVariant0.8.noRefCalls.vcf.gz
IN_FATHER_BACKGROUND_VCF=Polaris_VCF/ERR1955404_DeepVariant0.8.noRefCalls.vcf
FATHER_COVERAGE=45


# Unzip DeepVariant VCFs
gunzip -c $IN_PROBAND_BACKGROUND_VCFGZ > $IN_PROBAND_BACKGROUND_VCF
gunzip -c $IN_FATHER_BACKGROUND_VCFGZ > $IN_FATHER_BACKGROUND_VCF
gunzip -c $IN_MOTHER_BACKGROUND_VCFGZ > $IN_MOTHER_BACKGROUND_VCF

# Remove empty lines from VCF
sed -i '/^$/d' $MERGED_GB_VCF

# Split MERGED VCF to get sub-VCFs for each individual
bgzip -c $MERGED_GB_VCF > ${MERGED_GB_VCF}.gz
tabix ${MERGED_GB_VCF}.gz
file=${MERGED_GB_VCF}.gz
for sample in `bcftools query -l $file`; do
    bcftools view -c1 -Ov -s $sample -o ${file/.vcf*/.$sample.vcf} $file
done

# Modify variant format to match DeepVariant
# This will also output a background file saved as ${IN_PROBAND_BACKGROUND_VCF}_ReformatBackgound.vcf
#PROBAND
$PYTHON $REFORMAT_SCRIPT \
	-I $IN_PROBAND_GB_VCFA \
	-B $IN_PROBAND_BACKGROUND_VCF \
	-C $PROBAND_COVERAGE \
	-F proband \
	-R 

#MOTHER
$PYTHON $REFORMAT_SCRIPT \
        -I $IN_MOTHER_GB_VCFA \
        -B $IN_MOTHER_BACKGROUND_VCF \
        -C $MOTHER_COVERAGE \
        -F mother \
        -R

#FATHER
$PYTHON $REFORMAT_SCRIPT \
        -I $IN_FATHER_GB_VCFA \
        -B $IN_FATHER_BACKGROUND_VCF \
        -C $FATHER_COVERAGE \
        -F father \
        -R

# Organize your files
mv ${IN_PROBAND_BACKGROUND_VCF::-4}_Background_proband.vcf $WORKING_DIR/proband_Background.vcf
mv ${IN_PROBAND_GB_VCFA::-4}.vcf $WORKING_DIR/proband_PathoVar.vcf
mv ${IN_MOTHER_BACKGROUND_VCF::-4}_Background_mother.vcf $WORKING_DIR/mother_Background.vcf
mv ${IN_FATHER_BACKGROUND_VCF::-4}_Background_father.vcf $WORKING_DIR/father_Background.vcf

cd $WORKING_DIR
# BGZIP and tabix the new vcfs
bgzip -f proband_Background.vcf
tabix proband_Background.vcf.gz

bgzip -f proband_PathoVar.vcf
tabix proband_PathoVar.vcf.gz

bgzip -f mother_Background.vcf
tabix mother_Background.vcf.gz

bgzip -f father_Background.vcf
tabix father_Background.vcf.gz

# Concatenate VCFs for relevant individuals
bcftools concat -a proband_Background.vcf.gz proband_PathoVar.vcf.gz -O z -o proband_GeneBreaker.vcf.gz
tabix proband_GeneBreaker.vcf.gz

# Merge across all individuals
bcftools merge -0 proband_GeneBreaker.vcf.gz mother_Background.vcf.gz father_Background.vcf.gz -o Merged_GeneBreaker.vcf

bgzip Merged_GeneBreaker.vcf
tabix Merged_GeneBreaker.vcf.gz

