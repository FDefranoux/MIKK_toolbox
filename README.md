# Welcome to the MIKK panel Toolbox 

>
> Description: ...
>

## Installation

```
python3 -m venv venv_mikk_toolbox
source venv_mikk_toolbox/bin/activate
pip install -r requirement.txt
```

## Running the app

`streamlit run MIKK_toolbox.py`

## Parameters 

### General parameters

- [ ] Set up the initial session state.

#### Setting a working folder
- [ ] Using the `-f` option you can set up the folder containing the files you would like to work with. 

`streamlit run MIKK_toolbox.py -f finemapping`

- [ ] YAML to overwrite some of the file paths ?

#### YAML file 
- [ ] Using the `-y` option you can set up the parameters individually for the files.

Example of YAML file to setup the parameters of analysis 
>

#### Form ?

- [ ] ask to enter the file path manually ?

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
- [x] Implement general parameters (YAML and FOLDER + interactive and save to YAML)
- [x] Add ask for input if nothing exists 
- [ ] Check VEP idea
- [ ] Think how to link the different pages and their output together (selection of SNP → VEP analysis or send to FASTA etc)

