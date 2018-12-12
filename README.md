# Mendelian Disease Simulator Framework 

## Introduction 

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer elit urna, interdum ac augue in, efficitur consequat dui. Duis ac metus ut lorem molestie tincidunt. Nunc porttitor metus sit amet fringilla condimentum. Vestibulum sit amet ante mollis, viverra ligula ut, pharetra est. Suspendisse eu elit sit amet nibh volutpat ullamcorper. Vivamus ac lorem et massa dictum elementum. Duis tortor dui, tincidunt at urna at, venenatis hendrerit nunc. Etiam fringilla tristique odio, ut mattis mi. Sed dictum, enim vel interdum finibus, odio dui faucibus sapien, quis dignissim leo nulla in libero. Quisque eu iaculis eros. Phasellus malesuada mauris a elementum imperdiet. Nulla id felis in dolor auctor consectetur. Etiam nec tellus laoreet, ullamcorper mi accumsan, rhoncus velit.

## Directory Structure 

```
.
├── environment.yml
├── LICENCE
├── misc
├── README.md
├── sample_configs_vcfs
├── sample_ped
└── simulator
    ├── __init__.py
    ├── GUD
    ├── tests
    ├── scripts
    │   ├── cli.py
    │   ├── config_gen.py
    │   ├── family_vcf_gen.py
    │   └── __init__.py
    ├── clinvar.py
    ├── indel.py
    ├── short_tandem_repeat.py
    ├── single_nucleotide_variant.py
    ├── transcript.py
    ├── variant.py
    └── variants.py
```

## Dependancies 

- conda 
- linux/mac

## Modules 

1. Config Generator 

Under `simulator/scripts/config_gen.py` this module is used to to generate config files. This module can be called using the following command `python -m simulator.scripts.config_gen`. 

2. VCF Generator

Under `simulator/scripts/cli.py` this module is used to to generate vcf files from config files. This module can be called by the following command:

```
python -m simulator.scripts/cli -c <congif.json> -o <output>
## <config.json> is the config file, by default it will look for config.json .
## <output> is the name of the output, by default it will name the output: output.vcf .
```



3. Family Generator 

## Quickstart 

Follow these instructions to set up the package: 

```
git clone https://github.com/tamario/variant_simulation_framework.git
cd variant_simulation_framework
conda env create -f environment.yml
source activate var_sim
cd simulator/GUD
python setup.py install
cd ../../
python -m simulator.scripts.cli 
```

The output should be 3 files: output.vcf

## Running 

The simulator currently has two major points of entry, the configuration generator and the simulator. The standard workflow to run would be this 

1. Run configuration generator to generate config file that follows the EBNF. Instructions on what to enter are in the config_generator 

```
source activate var_sim
python -m simulator.scripts.config_gen
```

2. Run simulator to get vcfs

```
source activate var_sim
python -m simulator.scripts.cli -c <congig file name> -o <output file name>
```

## EBNF 

```
CONFIG ::= '{'GENE_INDEX','GENE_NAME', 'INHERITANCE', 'TRIO', 'VAR1', 'VAR2'}'

GENE_INDEX ::= number
GENE_NAME ::= inputstring
INHERITANCE ::= 'DE-NOVO' | 'BI-PARENTAL' | 'MATERNAL' | 'PATERNAL'
TRIO ::= 'TRIO' | 'SINGLE'
VAR1 ::= VARIANT_TEMPLATE
VAR2 ::= 'NONE' | VARIANT_TEMPLATE

VARIANT_TEMPLATE ::= '{'TYPE', 'REGION', 'LOCATION', 'IMPACT'}'

TYPE ::= 'SNV' | 'INDEL' | 'STR'
REGION ::= 'CODING' | 'UTR' | 'INTRONIC'
LOCATION ::= number | 'ANY' 
STR ::= {chrom, start, end} 
IMPACT ::=  
SNV_IMPACT_CODING := 'A' | 'T' | 'G' | 'C' | 'ANY' | ''MISSENSE | 'NONSENSE' | 'SILENT'
SNV_IMPACT_NONCODING := 'A' | 'T' | 'G' | 'C' | 'ANY'
INDEL_IMPACT ::= number
STR_IMPACT ::= number ## negative number is deletion of that size*motif size same logic for insertion 
```





