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

## Workflow

- Database creation
- Logic creation 
  1. select gene/phenotype + inheretance 
  2. select 1 or 2 variants
  3. use know n or design own variant 
  4. select var type (snv, indel, snv, str, sv, mei)
  5. select var impact (coding, utr, intronic, promoter, enhancer)
  6. present variants as vcf output 