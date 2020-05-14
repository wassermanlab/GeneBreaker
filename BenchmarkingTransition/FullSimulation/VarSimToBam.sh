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

# Where are you
SIMDIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/GeneBreaker/BenchmarkingTransition/FullSimulation/

# Activate environment
source $SIMDIR/../../opt/miniconda3/etc/profile.d/conda.sh
conda activate $SIMDIR/../../opt/GeneBreakerEnvironment

WORKING_DIR=$SIMDIR
RAW_DIR=$SIMDIR/GB_FullSimTest_out/

# concatenate the lanes together from varsim
cat $RAW_DIR/*.read1.fq.gz > $RAW_DIR/GB_FullSimTest.R1.fastq.gz
cat $RAW_DIR/*.read2.fq.gz > $RAW_DIR/GB_FullSimTest.R2.fastq.gz

FASTQR1=$RAW_DIR/GB_FullSimTest.R1.fastq.gz
FASTQR2=$RAW_DIR/GB_FullSimTest.R2.fastq.gz

ls $FASTQR1
ls $FASTQR2

#### GRCh37 ####
# map and convert against GRCh37
BWA_INDEX=/project/ex-ofornes-1/GENOME/GRCh37/GRCh37-lite.fa
SAMPLE_ID=GB_FullSimTest

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


