# Activate environment
source ../../opt/miniconda3/etc/profile.d/conda.sh
conda activate ../../opt/GeneBreakerEnvironment

# Set up this header 
# Software
PYTHON=python
REFORMAT_SCRIPT=$PWD/reformatSimToDeepVariant.py

# Run Details
GENE=WAS
GENOME=GRCh37
INHERITANCE=XLinkedRecessiveHemizygous
SEX=Male

# Make a new working directory for this simulation
WORKING_DIR=$PWD/${GENE}_${GENOME}_${INHERITANCE}_${SEX}
mkdir -p $WORKING_DIR

# Family-based settings
# Merged VCF file, here we pull from example directory
# If youre VCF is not named like this, make sure you change it to be.
# E.g. ../../Examples/MSH2_GRCh37_AutosomalDominantDeNovo.vcf 
MERGED_GB_VCF=$PWD/../../Examples/InheritanceTesting/${GENE}_${GENOME}_${INHERITANCE}_${SEX}.vcf
ls $MERGED_GB_VCF

if [ ! -f $MERGED_GB_VCF ];then
        echo "Variant VCF Not found"
        exit
fi

# Male proband; Family GB58 (AFR;GWD)
# Proband
IN_PROBAND_GB_VCF=${MERGED_GB_VCF::-4}.proband.vcf
IN_PROBAND_BACKGROUND_VCFGZ=$PWD/Polaris_VCF/$GENOME/ERR2304556_${GENOME}_DeepVariant_0.10.0_NoRefCalls.vcf.gz
IN_PROBAND_BACKGROUND_VCF=$PWD/Polaris_VCF/$GENOME/ERR2304556_${GENOME}_DeepVariant_0.10.0_NoRefCalls.vcf
PROBAND_COVERAGE=45

# Mother
IN_MOTHER_GB_VCF=${MERGED_GB_VCF::-4}.mother.vcf
IN_MOTHER_BACKGROUND_VCFGZ=$PWD/Polaris_VCF/$GENOME/ERR1955443_${GENOME}_DeepVariant_0.10.0_NoRefCalls.vcf.gz
IN_MOTHER_BACKGROUND_VCF=$PWD/Polaris_VCF/$GENOME/ERR1955443_${GENOME}_DeepVariant_0.10.0_NoRefCalls.vcf
MOTHER_COVERAGE=45

# Father
IN_FATHER_GB_VCF=${MERGED_GB_VCF::-4}.father.vcf
IN_FATHER_BACKGROUND_VCFGZ=$PWD/Polaris_VCF/$GENOME/ERR1955420_${GENOME}_DeepVariant_0.10.0_NoRefCalls.vcf.gz
IN_FATHER_BACKGROUND_VCF=$PWD/Polaris_VCF/$GENOME/ERR1955420_${GENOME}_DeepVariant_0.10.0_NoRefCalls.vcf
FATHER_COVERAGE=45


# Unzip DeepVariant VCFs
gunzip -f -c $IN_PROBAND_BACKGROUND_VCFGZ > $IN_PROBAND_BACKGROUND_VCF
gunzip -f -c $IN_FATHER_BACKGROUND_VCFGZ > $IN_FATHER_BACKGROUND_VCF
gunzip -f -c $IN_MOTHER_BACKGROUND_VCFGZ > $IN_MOTHER_BACKGROUND_VCF

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
	-I $IN_PROBAND_GB_VCF \
	-B $IN_PROBAND_BACKGROUND_VCF \
	-C $PROBAND_COVERAGE \
	-F proband 

#MOTHER
$PYTHON $REFORMAT_SCRIPT \
        -I $IN_MOTHER_GB_VCF \
        -B $IN_MOTHER_BACKGROUND_VCF \
        -C $MOTHER_COVERAGE \
        -F mother 

#FATHER
$PYTHON $REFORMAT_SCRIPT \
        -I $IN_FATHER_GB_VCF \
        -B $IN_FATHER_BACKGROUND_VCF \
        -C $FATHER_COVERAGE \
        -F father 

# Organize your files
mv ${IN_PROBAND_BACKGROUND_VCF::-4}_Background_proband.vcf $WORKING_DIR/proband_Background.vcf
mv ${IN_PROBAND_GB_VCF::-4}.vcf $WORKING_DIR/proband_PathoVar.vcf
mv ${IN_MOTHER_BACKGROUND_VCF::-4}_Background_mother.vcf $WORKING_DIR/mother_Background.vcf
mv ${IN_MOTHER_GB_VCF::-4}.vcf $WORKING_DIR/mother_PathoVar.vcf
mv ${IN_FATHER_BACKGROUND_VCF::-4}_Background_father.vcf $WORKING_DIR/father_Background.vcf
mv ${IN_FATHER_GB_VCF::-4}.vcf $WORKING_DIR/father_PathoVar.vcf

cd $WORKING_DIR
# BGZIP and tabix the new vcfs
bgzip -f proband_Background.vcf
tabix proband_Background.vcf.gz

bgzip -f proband_PathoVar.vcf
tabix proband_PathoVar.vcf.gz

bgzip -f mother_Background.vcf
tabix mother_Background.vcf.gz

bgzip -f mother_PathoVar.vcf
tabix mother_PathoVar.vcf.gz

bgzip -f father_Background.vcf
tabix father_Background.vcf.gz

bgzip -f father_PathoVar.vcf
tabix father_PathoVar.vcf.gz

# Concatenate VCFs for relevant individuals
bcftools concat -a proband_Background.vcf.gz proband_PathoVar.vcf.gz -O z -o proband_GeneBreaker.vcf.gz
tabix proband_GeneBreaker.vcf.gz

bcftools concat -a mother_Background.vcf.gz mother_PathoVar.vcf.gz -O z -o mother_GeneBreaker.vcf.gz
tabix mother_GeneBreaker.vcf.gz

bcftools concat -a father_Background.vcf.gz father_PathoVar.vcf.gz -O z -o father_GeneBreaker.vcf.gz
tabix father_GeneBreaker.vcf.gz

# Merge across all individuals
bcftools merge -0 proband_GeneBreaker.vcf.gz mother_GeneBreaker.vcf.gz father_GeneBreaker.vcf.gz -o Merged_GeneBreaker.vcf
bgzip -f Merged_GeneBreaker.vcf
tabix Merged_GeneBreaker.vcf.gz


# Finished!
echo "Finished!"
