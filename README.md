# ARGfinder
A tool for for large-scale analysis of antimicrobial resistance genes and their flanking regions in metagenomic datasets.


## Introduction

ARGfinder is a newly developed snakemake pipeline designed to analyze read distances, abundancies and genomic flanking regions of ARGs in metagenomic sequencing data. It has been adapted to work for short-read sequencing datasets and includes the recently made Panres database, a combined collection of current ARG databases, and ARGextender, an assembly tool for extending the genomic flanking region arouns genes of interest.

ARGfinder uses the following tools:


* ``` enaBrowserTools ``` for downloading raw reads from ENA
* ``` fastp ``` for trimming and QC of raw reads
* ``` KMA ``` for alignment of raw reads against reference databases
* ``` ARGextender ``` for extracting the genomic flanking regions around ARGs
* ``` Mash ``` for creating sketches to estimate genetic distances


The tool is described in *Paper_here*

## Getting started

The best way to install the ARGfinder pipeline is to clone this github repository. The pipeline uses the Conda package manager to deploy the defined software packages in the specified version without requiring any admin/root priviledges.

```
git clone https://github.com/genomicepidemiology/ARGfinder.git
```
This command will create the directory ARGfinder in the current directory.

Since ARGfinder is a snakemake pipeline the user should have snakemake installed. 