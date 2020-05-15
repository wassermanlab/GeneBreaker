# Optional Header info for simulating using PBS
#!/bin/bash
#PBS -A ex-ofornes-1
#PBS -V
#PBS -N Exomiser
#PBS -m bea
#PBS -M prichmond@cmmt.ubc.ca
#PBS -l select=1:ncpus=8:mem=30gb
#PBS -l walltime=4:0:0

NSLOTS=8

# get to proper directory
cd /scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/GeneBreaker/BenchmarkingTransition/BuryVariant/ExomiserTesting/

#ABCD1
EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=ABCD1_GRCh38_XLinkedHemizygousRecessiveDeNovo_Male.yml

java -jar $EXOMISER_DIR/$EXOMISER \
	 --analysis $CONFIG_YML \
	--spring.config.location=$EXOMISER_DIR

#CFTR
EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=CFTR_GRCh38_AutosomalRecessiveCompoundHeterozygous_Female.yml

java -jar $EXOMISER_DIR/$EXOMISER \
	 --analysis $CONFIG_YML \
	--spring.config.location=$EXOMISER_DIR

#INPP5E
EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=INPP5E_GRCh38_AutosomalRecessiveCompoundHeterozygousDeNovo_Female.yml

java -jar $EXOMISER_DIR/$EXOMISER \
	 --analysis $CONFIG_YML \
	--spring.config.location=$EXOMISER_DIR

#JAK1
EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=JAK1_GRCh37_AutosomalDominant_Male.yml

java -jar $EXOMISER_DIR/$EXOMISER \
	 --analysis $CONFIG_YML \
	--spring.config.location=$EXOMISER_DIR

#MALT1
EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=MALT1_GRCh37_AutosomalRecessiveHomozygous_Female.yml

java -jar $EXOMISER_DIR/$EXOMISER \
	 --analysis $CONFIG_YML \
	--spring.config.location=$EXOMISER_DIR

#MECP2
EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=MECP2_GRCh38_XLinkedDominantDeNovo_Female.yml

java -jar $EXOMISER_DIR/$EXOMISER \
	 --analysis $CONFIG_YML \
	--spring.config.location=$EXOMISER_DIR
#MSH2
EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=MSH2_GRCh37_AutosomalDominantDeNovo_Male.yml

java -jar $EXOMISER_DIR/$EXOMISER \
	 --analysis $CONFIG_YML \
	--spring.config.location=$EXOMISER_DIR

#SLC6A8
EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=SLC6A8_GRCh38_XLinkedRecessiveCompoundHeterozygousDeNovo_Female.yml

java -jar $EXOMISER_DIR/$EXOMISER \
	 --analysis $CONFIG_YML \
	--spring.config.location=$EXOMISER_DIR

#SRY
EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=SRY_GRCh38_YLinkedDominantDeNovo_Male.yml

java -jar $EXOMISER_DIR/$EXOMISER \
	 --analysis $CONFIG_YML \
	--spring.config.location=$EXOMISER_DIR

#WAS
EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=WAS_GRCh37_XLinkedRecessiveHemizygous_Male.yml

java -jar $EXOMISER_DIR/$EXOMISER \
	 --analysis $CONFIG_YML \
	--spring.config.location=$EXOMISER_DIR
