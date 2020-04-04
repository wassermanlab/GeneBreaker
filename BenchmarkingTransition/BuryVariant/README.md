# Bury Variants

> The purpose of this pipeline is to bury SNV/indel variants within VCF files

## Overview

1. Get necessary tools by deploying a miniconda within the main repository directory, which will end up under 'opt'

```
sh Install_Downstream_Tools.sh
```

2. Download DeepVariant VCFs from Zenodo

```
sh GetBackgroundVars.sh
```

3. Run an example, MSH2.

```
sh IntegrateVariants_MSH2_GRCh37_AutosomalDominantDeNovo.sh
```

4. (optional) Run your own test cases!


## Building your own embeddings

This is a stepwise process. Assume for this example, that I've just simulated a de novo variant in the gene PKD1 in GRCh38.

1. Select the inheritance mode script to match your simulated inheritance mode. Several examples exist here.

```
cp IntegrateVariants_AutosomalDominantDeNovo.sh PKD1_GRCh38_AutosomalDominantDeNovo.sh
```

2. Edit the script 'IntegrateVariants\_PKD1\_GRCh38\_AutosomalDominantDeNovo.sh' Change the filepaths, and genome version, as directed in the shell script to match up to your input simulated VCF. If you're a boss then you'll use vim.

```
vim IntegrateVariants_PKD1_GRCh38_AutosomalDominantDeNovo.sh
```

3. Execute the variant burying.

```
sh  IntegrateVariants_PKD1_GRCh38_AutosomalDominantDeNovo.sh 
```

## Optional changes:

+ If you prefer to use your own VCFs for embedding, then you will have to make sure that the variant line and formatting matches between your simulated VCF variantand the VCFs you wish to use for embedding. We take care of this within a reformatSimToDeepVariant.py script, which could be replicated for your format of choice.




