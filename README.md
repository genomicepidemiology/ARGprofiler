# ARGfinder
A tool for for large-scale analysis of antimicrobial resistance genes and their flanking regions in metagenomic datasets.

<img src="ARGfinder_pipeline.png" alt="ARGfinder pipeline">

## Introduction

ARGfinder is a newly developed Snakemake pipeline designed to analyze read distances, abundancies and genomic flanking regions of ARGs in metagenomic sequencing data. It has been adapted to work for short-read sequencing datasets and includes the recently made Panres database, a combined collection of current ARG databases, and ARGextender, an assembly tool for extending the genomic flanking region arouns genes of interest.

ARGfinder uses the following tools:


* ``` enaBrowserTools ``` for downloading raw reads from ENA
* ``` fastp ``` for trimming and QC of raw reads
* ``` KMA ``` for alignment of raw reads against reference databases
* ``` ARGextender ``` for extracting the genomic flanking regions around ARGs
* ``` Mash ``` for creating sketches to estimate genetic distances


The tool is described in *Paper_here*

## Installation

The best way to install the ARGfinder pipeline is to clone this github repository. The pipeline uses the Conda package manager to deploy the defined software packages in the specified version without requiring any admin/root priviledges.

```
git clone https://github.com/genomicepidemiology/ARGfinder.git
```
This command will create the directory ARGfinder in the current directory.

Since ARGfinder is a Snakemake pipeline the user should have Snakemake workflow management installed. 

## Getting started

There are some prerequisites for using ARGfinder:

* The user needs to download the 2 reference databases (mOTUs and PanRes), place them in the the correct directory and then index them with KMA. 

	* Create mOTUs and panres subdirectories in the ``` prerequisites ``` directory.

	* For ``` mOTUs``` the user has to download the database from zenodo, unzip it and then index it with KMA:
		1. Move to mOTUs directory in the ``` prerequisites ``` directory
		2. wget https://zenodo.org/record/5140350/files/db_mOTU_v3.0.1.tar.gz
		3. tar -xzf db_mOTU_v3.0.1.tar.gz
		4. kma index -i mOTUs/db_mOTU_DB_CEN.fasta -o db_mOTU (For KMA instructions you can check  <a href="https://bitbucket.org/genomicepidemiology/kma/src/master/">KMA</a>)

	* For ``` PanRes``` the user has to download the database from zenodo, unzip it and then index it with KMA:
		1. Move to panres directory in the ``` prerequisites ``` directory
		2. wget https://zenodo.org/record/
		3. tar -xzf 
		4. kma index -i panres/panres_db -o db_panres (For KMA instructions you can check  <a href="https://bitbucket.org/genomicepidemiology/kma/src/master/">KMA</a>)

* The pipeline make use of Snakemake profiles to specify configuration of the pipeline. The various flags that are needed, are specified in the files of the ``` profile_argfinder ``` directory.

## Input

ARGfinder takes as input a file with the following format:
```
{"ERR3593315":{"type":"PAIRED"},"SRR7533096":{"type":"SINGLE"}}
```

The user has aldo to specify the names of the file in the Snakefile (with open...)

## Running ARGfinder

## Output