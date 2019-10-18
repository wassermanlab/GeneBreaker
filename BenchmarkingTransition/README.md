# Benchmarking Transition

> This directory is for the addition of downstream benchmarking pipelines which take as input the VCF(s) from the simulation.

## Bury Variant Within DeepVariant

> This approach is designed for testing interpretation of SNV/indels

1. Get necessary tools for adding variant (VCFTools). 
	
	$ sh GetTools.sh

2. Collect DeepVariant VCFs (need link for this, maybe host them on Zenodoo)

	$ sh GetVariantVCFs.sh

3. Use BuryVCF.sh

	$ sh BuryVCF.sh

## Simulate Variant using VarSim

> This approach is designed to facilitate the simulation of 




