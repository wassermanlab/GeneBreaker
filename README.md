# Mendelian Disease Simulator Framework 

## Introduction 



## Directory Structure 

```
./variant_simulation_framework
├── README.md
├── requirements.txt
├── simulator
│   ├── gene.py     #class
│   ├── variants.py #class
│   └── variant_types #modules 
│       ├── copy_number_variant.py
│       ├── indel.py
│       ├── mobile_element.py
│       ├── single_nucleotide_variant.py 
│       ├── single_tandem_repeat.py
│       └── structural_variant.py
├── tests
└── varsim
    ├── run_varsim.sh
    └── setup.sh

```

## EBNF 

NEW EBNF

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



## Notes/References

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4422378/table/Tab2/?report=objectonly











