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

```

## gene
GENE ::='{'SEQ', 'CHR', 'START_REG', 'STOP_REG', 'ELEMENTS'}'
SEQ ::= string
CHR ::= number(1:22) | X | Y
START_REG  ::= number
STOP_REG  ::= number
ELEMENTS  ::= 'ELEMENTS:[' (ELEMENT ',')* ELEMENT ']'
ELEMENT ::= [EL_TYPE, START_POS, STOP_POS]
EL_TYPE ::= 'EXON' | 'UTR' | 'INTRON' | 'PROMOTER' | 'ENHANCER
START_POS ::= number
STop_POS ::= number

#variants 
VARIANTS ::= '{'DESCRIPTION', 'INHERITANCE', 'VAR1', 'VAR2'}'
INHERITANCE ::= 'DE-NOVO' | 'BI-PARENTAL' | 'MATERNAL' | 'PATERNAL'
DESCRIPTION ::= optional string discription 
VAR1 ::= 'NONE' | VARIANT_TEMPLATE
VAR2 ::= 'NONE' | VARIANT_TEMPLATE

## Variant template
VARIANT_TEMPLATE ::='{'TYPE', 'REGION', 'IMPACT'}'

TYPE ::= 'TYPE: 'SNV' | 'INDEL' | 'CNV' | 'SV' | 'STR' | 'MEI''
REGION ::= 'REGION: CODING' | 'UTR' | 'INTRONIC' | 'PROMOTER' | 'ENHANCER''
IMPACT ::= 'IMPACT: {'TYPE_IMPACT', 'LOCATION'}'
TRIO ::= 'True' | 'False'

SNV_IMPACT ::= 'MISSENSE' | 'NONSENSE' | 'SILENT'
INDEL_IMPACT ::= NUMBER ## represents the number of bases to delete or insert where positive numbers are insertions and negative deletions
CNV_IMPACT ::= 
SV_IMPACT ::= 
STR_IMPACT ::= 
MEI_IMPACT ::= 

LOCATION ::= 'ANY' | number ## where number is position on chromosome 

## examples 
{'TYPE': 'SNV',
'REGION' 'CODING',
'IMPACT': {'SNV_IMPACT': 'MISSENSE',
			'LOCATION': ANY}}

{'TYPE': 'INDEL',
'REGION' 'INTRONIC',
'IMPACT': {'INDEL_IMPACT': -50,
			'LOCATION': 11756}}
```





## Workflow

- Database creation
- Logic creation 
  1. select gene/phenotype + inheretance 
  2. select 1 or 2 variants
  3. use know n or design own variant 
  4. select var type (snv, indel, snv, str, sv, mei)
  5. select var impact (coding, utr, intronic, promoter, enhancer)
  6. present variants as vcf output 

## Notes/References

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4422378/table/Tab2/?report=objectonly











