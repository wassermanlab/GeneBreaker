# Bury Variants

> The purpose of this pipeline is to bury SNV/indel variants within VCF files

## Overview

1. Get necessary tools by deploying a miniconda within the *main repository directory*, which will end up under 'GeneBreaker/opt/':

```
sh Install_Downstream_Tools.sh
```

2. Download DeepVariant VCFs from Zenodo, and unwind them into the Polaris\_VCF/GRCh37/ and Polaris\_VCF/GRCh38/ directories. This is handled with a single script:

```
sh GetBackgroundVars.sh
```

3. Run an example, IntegrateVariants\_Male\_GRCh37.sh.

```
sh IntegrateVariants_Male_GRCh37.sh
```  

4. (optional) Run your own test cases! (You'll build your own below)  


## Building your own case

This is a stepwise process. Assume for this example, that I've just simulated an autosomal recessive compound heterozygous ~de novo~ variant in SLC6A8, within the GRCh38 reference genome.
I'll copy the output VCF and output PED from GeneBreaker to be formatted like this:
```
SLC6A8_GRCh38_XLinkedRecessiveCompoundHeterozygousDeNovo.vcf
SLC6A8_GRCh38_XLinkedRecessiveCompoundHeterozygousDeNovo.ped
```

1. Select the matching sex and genome integration script. 

```
cp IntegrateVariants_ GRCh38 PKD1\_GRCh38\_AutosomalDominantDeNovo.sh
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


