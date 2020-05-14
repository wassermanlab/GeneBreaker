# Phillip Richmond
# Optional Header info for simulating using PBS
#!/bin/bash
#PBS -A ex-ofornes-1
#PBS -V
#PBS -N VarSimToBam
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
#PBS -l select=1:ncpus=32:mem=150gb
#PBS -l walltime=48:0:0

NSLOTS=30

# Activate environment
SIMDIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/GeneBreaker/BenchmarkingTransition/FullSimulation/
source $SIMDIR/varsim/opt/miniconda2/etc/profile.d/conda.sh
#conda activate $SIMDIR/varsim/opt/GeneBreakerEnvironment
conda activate $SIMDIR/varsim/opt/miniconda2/
cd $SIMDIR

# Varsim variables
VARSIM_DIR=$PWD/varsim/
OPT_DIR=${VARSIM_DIR}/opt
SAMPLE_ID=GB_LargeVarsMEIsDarkCamo

# Test run varsim to generate 1x coverage data
#export PATH=${OPT_DIR}/jdk1.8.0_131/bin:$PATH
#$SIMDIR/varsim/varsim.py --vcfs $SIMDIR/../../Examples/InheritanceTesting/JAK1_GRCh37_AutosomalDominant_Male.vcf \
#	$SIMDIR/../../Examples/LargeVariants/DMD_GRCh37_XLinkedRecessiveCompoundHeterozygousCNV_Female.vcf \
#	$SIMDIR/../../Examples/LargeVariants/GLS_GRCh37_AutosomalRecessiveHomozygousSTR_Female.vcf \
#	$SIMDIR/../../Examples/LargeVariants/SLC2A1_GRCh37_AutosomalDominantDeNovoCNV_Female.vcf \
#	$SIMDIR/../../Examples/LargeVariants/IKBKG_GRCh37_XLinkedDominantDeNovoMEI_Female_forVarSim.vcf \
#	$SIMDIR/../../Examples/LargeVariants/SLC17A5_GRCh37_AutosomalRecessiveHomozygousMEI_Female_forVarSim.vcf \
#	$SIMDIR/../../Examples/DarkCamo/CFC1_GRCh37_AutosomalDominantDeNovoDarkCamo_Female.vcf \
#	$SIMDIR/../../Examples/DarkCamo/RPGR_GRCh37_XLinkedRecessiveCompoundHeterozygousDarkCamo_Female_forVarSim.vcf \
#	$SIMDIR/../../Examples/DarkCamo/MAF_GRCh37_AutosomalDominantDeNovoDarkCamo_Female.vcf \
#	$SIMDIR/../../Examples/DarkCamo/SMN2_GRCh37_AutosomalRecessiveHomozygousDarkCamo_Female.vcf \
#	--vc_in_vcf 00-common_all.vcf \
#	--sv_insert_seq insert_seq.txt \
#	--sv_dgv GRCh37_hg19_supportingvariants_2013-07-23.txt \
#	--reference hs37d5.fa --id simu --read_length 150 --vc_num_snp 3000000 --vc_num_ins 100000 \
#	--vc_num_del 100000 --vc_num_mnp 5 --vc_num_complex 5 --sv_num_ins 2000 \
#	--sv_num_del 2000 --sv_num_dup 200 --sv_num_inv 1000 --sv_percent_novel 0.01 \
#	--vc_percent_novel 0.01 --mean_fragment_size 450 --sd_fragment_size 100 \
#	--vc_min_length_lim 0 --vc_max_length_lim 49 --sv_min_length_lim 50 \
#	--sv_max_length_lim 1000000 --nlanes 3 --total_coverage 30 \
#	--java_max_mem 50g \
#	--simulator_executable ${OPT_DIR}/ART/art_bin_VanillaIceCream/art_illumina --out_dir ${SAMPLE_ID}_out \
#	--log_dir ${SAMPLE_ID}_log --work_dir ${SAMPLE_ID}_work \
#	--simulator art 
#
#exit
## Part II: Map & Convert

# Where are you
SIMDIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/GeneBreaker/BenchmarkingTransition/FullSimulation/

# Activate environment
source $SIMDIR/../../opt/miniconda3/etc/profile.d/conda.sh
conda activate $SIMDIR/../../opt/GeneBreakerEnvironment

WORKING_DIR=$SIMDIR
RAW_DIR=$SIMDIR/${SAMPLE_ID}_out/

# concatenate the lanes together from varsim
cat $RAW_DIR/*.read1.fq.gz > $RAW_DIR/GB_LargeVarsMEIsDarkCamo.R1.fastq.gz
cat $RAW_DIR/*.read2.fq.gz > $RAW_DIR/GB_LargeVarsMEIsDarkCamo.R2.fastq.gz

FASTQR1=$RAW_DIR/GB_LargeVarsMEIsDarkCamo.R1.fastq.gz
FASTQR2=$RAW_DIR/GB_LargeVarsMEIsDarkCamo.R2.fastq.gz

ls $FASTQR1
ls $FASTQR2

#### GRCh37 ####
# map and convert against GRCh37
BWA_INDEX=/project/ex-ofornes-1/GENOME/GRCh37/GRCh37-lite.fa

cd $WORKING_DIR

# Map reads
bwa mem $BWA_INDEX \
        -t $NSLOTS \
        -R "@RG\tID:$SAMPLE_ID\tSM:$SAMPLE_ID\tPL:illumina" \
        -M \
        $FASTQR1 \
        $FASTQR2 \
        > $WORKING_DIR$SAMPLE_ID.sam

# convert to binary and index
samtools view -@ $NSLOTS -ubS $WORKING_DIR$SAMPLE_ID'.sam' \
        | samtools sort - -@ $NSLOTS  -T $WORKING_DIR$SAMPLE_ID'.sorted' -O BAM -o $WORKING_DIR$SAMPLE_ID'.sorted.bam'
samtools index $WORKING_DIR$SAMPLE_ID'.sorted.bam'

rm $WORKING_DIR$SAMPLE_ID'.sam'

