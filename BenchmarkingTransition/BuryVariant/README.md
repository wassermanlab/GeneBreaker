# Bury Variants

> The purpose of this pipeline is to bury SNV/indel variants within VCF files

## Overview

1. Get necessary tools by deploying a miniconda within this directory, which will end up under 'opt'

	sh GetTools.sh

2. Download DeepVariant VCFs from Zenodo




## Optional changes:

+ If you prefer to use your own VCFs for embedding, then you will have to make sure that the variant line and formatting matches between your simulated VCF variantand the VCFs you wish to use for embedding. We take care of this within a reformatSimToDeepVariant.py script, which could be replicated for your format of choice.




