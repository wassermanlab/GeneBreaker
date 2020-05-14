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

EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=JAK1_GRCh37_AutosomalDominant_Male.yml

java -jar $EXOMISER_DIR/$EXOMISER \
	 --analysis $CONFIG_YML \
	--spring.config.location=$EXOMISER_DIR
