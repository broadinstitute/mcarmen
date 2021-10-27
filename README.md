# mCARMEN
This repo contains relevant codebase for CARMEN-RVP and CARMEN-VARIANT analyses.

# CARMEN-RVP analysis tool
CARMEN-RVP is a CRISPR/Cas13 based assay detecting 9 respiratory viruses which is run on a Fluidigm Biomark HD instrument. The readout is raw fluorescence data extracted from the Fluidigm RT-PCR software. The exported data can then be further analyzed with this tool.

src folder:\
The main script for the analysis is called 'carmen_rvp_analysis.py'. The other scripts contain functions that will be imported during the analysis run.

img folder:\
Icons for graphical user interface.

archive-expertversion folder:\
For R&D use only. This contains an older version of the tool that provides the user with more flexibility on the input.

## Requirements
Python dependencies:
- pandas==1.0.5
- seaborn==0.10.1
- matplotlib==3.2.2
- Gooey==1.0.7
- numpy==1.18.5
- xlrd==1.1.0

Inputs:
- Fluidigm Biomark raw data file (.csv)
- Sample assignment sheet (.xlsx): Template "assignment-template.xlsx" can be downloaded here.

## Installation

Clone this github repository to always have the most up-to-date version. An executable version for Windows and Mac will be coming soon.

## Preparation of the assignment sheet
The assignment sheet can be downloaded from this repository. \
Samples controls have to be named 'EC', 'NTC', 'NDC', 'CPC'. \
The control assay without crRNA has to be called 'no-crRNA' and the assay for RNaseP has to be called 'RNaseP' or 'RNAseP'.

## Workflow

1. If the command line is used, run the the main script with the following command:\
$ pythonw ./src/carmen_rvp_analysis.py
2. Choose an output prefix and an output directory
3. Upload the assignment sheet and the raw data file
4. Click start run.

# CARMEN-VARIANT analysis tool
CARMEN-VARIANT is a CRISPR/Cas13 based assay detecting COVID-19 variants and SNP mutations of concern which is run on a Fluidigm Biomark HD instrument. The readout is raw fluorescence data extracted from the Fluidigm RT-PCR software. The exported data can then be further analyzed with this tool.

src folder:\
The main scripts for the analysis:
- carmen_variant_dataparser.py: will parse the data exported from the Fluidigm RT-PCR software
- carmen_variant_analysis.py: will determine the test positivity for each SNP and the ultimate variant call

## Requirements
Python dependencies:
- pandas==1.0.5
- seaborn==0.10.1
- matplotlib==3.2.2
- numpy==1.18.5

Inputs:
- Fluidigm Biomark raw data file (.csv)
- Sample assignment sheet (.xlsx): Template "assignment-template.xlsx" can be downloaded here.
- FCthresholds (.csv): contains the thresholds for each crRNA guide pairs for a single SNP.

## Installation

Clone this github repository to always have the most up-to-date version. An executable version for Windows and Mac will be coming soon.

## Preparation of the assignment sheet
The assignment sheet can be downloaded from this repository. \
Samples controls have to be named 'EC', 'NTC', 'NDC', 'CPC'. \
The control assay without crRNA has to be called 'no-crRNA' and the assay for RNaseP has to be called 'RNaseP' or 'RNAseP'.

## Workflow

$ python ./src/carmen_variant_dataparser.py *arglist
$ python ./src/carmen_variant_analysis.py *arglist
