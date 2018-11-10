# Mendelian Disease Simulator Framework 

## Introduction 

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer elit urna, interdum ac augue in, efficitur consequat dui. Duis ac metus ut lorem molestie tincidunt. Nunc porttitor metus sit amet fringilla condimentum. Vestibulum sit amet ante mollis, viverra ligula ut, pharetra est. Suspendisse eu elit sit amet nibh volutpat ullamcorper. Vivamus ac lorem et massa dictum elementum. Duis tortor dui, tincidunt at urna at, venenatis hendrerit nunc. Etiam fringilla tristique odio, ut mattis mi. Sed dictum, enim vel interdum finibus, odio dui faucibus sapien, quis dignissim leo nulla in libero. Quisque eu iaculis eros. Phasellus malesuada mauris a elementum imperdiet. Nulla id felis in dolor auctor consectetur. Etiam nec tellus laoreet, ullamcorper mi accumsan, rhoncus velit.

## Directory Structure 

```
./variant_simulation_framework
├── config.json
├── environment.yml
├── LICENCE
├── output.vcf
├── README.md
├── requirements.txt
├── simulator
│   ├── GUD ...
│   ├── tests ...
│   ├── scripts
│   │   ├── __init__.py
│   │   ├── config_gen.py
│   │   └── cli.py
│   ├── __init__.py
│   ├── indel.py
│   ├── single_nucleotide_variant.py
│   ├── transcript.py
│   ├── variant.py
│   └── variants.py
└── varsim
    ├── run_varsim.sh
    └── setup.sh


```

## Dependancies 

- conda 

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

The output should be 3 files: ouput.vcf.child, ouput.vcf.mother, and ouput.vcf.father

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
# CONFIG
CONFIG ::= '{'GENE_INDEX','GENE_NAME', 'INHERITANCE', 'TRIO', 'VAR1', 'VAR2',}'

GENE_INDEX ::= INT
GENE_NAME ::= STRING
INHERITANCE ::= 'DE-NOVO' | 'BI-PARENTAL' | 'MATERNAL' | 'PATERNAL'
TRIO ::= 'TRIO' | 'SINGLE'
VAR1 ::= VARIANT_TEMPLATE
VAR2 ::= 'NONE' | VARIANT_TEMPLATE

# VARIANT_TEMPLATE
VARIANT_TEMPLATE ::= '{'TYPE', 'REGION', 'LOCATION', 'IMPACT'}'

TYPE ::= 'SNV' | 'INDEL' 
REGION ::= 'CODING' | 'UTR' | 'INTRONIC'
LOCATION ::= INT | 'ANY'
IMPACT ::= multiple 
	- SNV_IMPACT_CODING := 'A' | 'T' | 'G' | 'C' | 'ANY' | ''MISSENSE | 'NONSENSE' | 'SILENT'
	- SNV_IMPACT_NONCODING := 'A' | 'T' | 'G' | 'C' | 'ANY'
	- INDEL_IMPACT ::= INT
```





