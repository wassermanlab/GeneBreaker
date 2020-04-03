EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/scratch/st-wasserww-1/GeneBreaker/BenchmarkingTransition/BuryVariant/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=/scratch/st-wasserww-1/GeneBreaker/BenchmarkingTransition/BuryVariant/GeneBreaker_DeNovo_genome.yml 

java -jar $EXOMISER_DIR/$EXOMISER \
	 --analysis $CONFIG_YML \
	--spring.config.location=/scratch/st-wasserww-1/GeneBreaker/BenchmarkingTransition/BuryVariant/Exomiser/exomiser-cli-12.1.0/
