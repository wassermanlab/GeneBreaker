EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/scratch/ex-ofornes-1/RICHMOND/GENEBREAKER/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=MSH2_GRCh37_AutosomalDominantDeNovo.yml

java -jar $EXOMISER_DIR/$EXOMISER \
	 --analysis $CONFIG_YML \
	--spring.config.location=$EXOMISER_DIR
