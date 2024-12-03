# Readme

>
> Description: This scripts are made to have a little GUI to explore data adapted for our MIKK panel but should work for similar dataset in other organisms.
> 1. FASTA retriever: Convert VCF files in alignement to a reference. You will need the fasta files per chr of your reference.
> 2. GWAS analysis: mostly adapted to analyse dataset after runing the flexlmm pipeline
> 3. Region information: Small visualisation for VEP files
>

## Installation

```
python3 -m venv venv_mikk_toolbox
source venv_mikk_toolbox/bin/activate
pip install -r requirement.txt
```

## Running the app

`streamlit run MIKK_toolbox.py`

`streamlit run MIKK_toolbox.py -- -y path/to/yaml/params.yml`

## Parameters 

### General parameters

-y : YAML parameters path. Not mandatory can be added through the app.


### Function description and their parameters 

#### FASTA retriever

| variable name   | required    | default value |
|-----------------|-------------|---------------|
| vcf_file        | True        | None          |
| Fasta_directory | True        | None          |
| Cram file       | False       | None          |
| chr             | Interactive | None          |
| start           | Interactive | None          |
| end             | Interactive | None          |


#### GWAS analysis

> At the moment this script is very specific to the use of flexlmm!

| variable name                | required              | default value |
|------------------------------|-----------------------|---------------|
| VCF file                     | True                  | None          |
| **flexlmm output directory** | True                  | None          |
| phenotype file               | **Only in flexlmm ?** | None          |
| covariate file               | **Only in flexlmm ?** | None          |
| chr                          | Interactive           | None          |


#### Region information

| variable name           | required    | default value |
|-------------------------|-------------|---------------|
| info_file (format gff3) | True        | None          |
| chr                     | Interactive | None          |
| start                   | Interactive | None          |
| end                     | Interactive | None          |
| marker pos              | False       | None          |



#### TODOs
- [ ] Extract parameters directly from flexlmm param file
- [ ] Think how to link the different pages and their output together (selection of SNP → VEP analysis or send to FASTA etc)


