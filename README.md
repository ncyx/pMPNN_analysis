# Automated Analysis of ProteinMPNN Designed Sequences (pMPNN_analysis)
![pipeline](/figure_pipeline/pipeline.png)
## Description

This repository contains scripts on automated data analysis and visualization of biophysical parameters of pMPNN-designed (benchmarked for monomeric) proteins. 

The main usage was meant to analyze biophysical parameters that likely contribute to the enhanced protein solubility: (i) hydrophilicity and (ii) charge state of designed solvent-exposed residues, and (iii) *in silico* formed core size. Developed scripts extract from the output .fa and predicted structures of designed variants (as .pdb) plddt, rmsd to native, and the above-mentioned biophysical properties.

* [examples](https://github.com/ncyx/pMPNN_analysis/tree/main/examples) contains example scripts and raw data (fasta and pdb) used to design proteins with pMPNN and analyze the results.

## Installation 
Installation can be carried out directly from the .yml:

`conda env create -f env_dependencies.yml`

## Usage 

* In the [scripts](https://github.com/ncyx/pMPNN_analysis/tree/main/scripts) folder, there is a JupyterNotebook that illustrates the pipeline for data analysis.
* To go through this tutorial, you should clone the repository and make sure that the paths match the location of scripts:

  `git clone https://github.com/ncyx/pMPNN_analysis.git`
  
  `cd pmPNN_analysis/scripts`
