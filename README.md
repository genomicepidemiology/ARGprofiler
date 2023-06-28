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
		4. kma index -i db_mOTU/db_mOTU_DB_CEN.fasta -o db_mOTU (For KMA instructions you can check  <a href="https://bitbucket.org/genomicepidemiology/kma/src/master/">KMA</a>)

	* For ``` PanRes``` the user has to download the database from zenodo, unzip it and then index it with KMA:
		1. Move to panres directory in the ``` prerequisites ``` directory
		2. wget https://zenodo.org/record/
		3. tar -xzf 
		4. kma index -i pan.fa -o panres_db (For KMA instructions you can check  <a href="https://bitbucket.org/genomicepidemiology/kma/src/master/">KMA</a>)

* The pipeline make use of Snakemake profiles to specify configuration of the pipeline. The various flags that are needed, are specified in the files of the ``` profile_argfinder ``` directory.
	
	* We provide a config file for executing the pipeline in an HPC with qsub

## Input

ARGfinder takes as input a file with the following format:

```
{Accession_run_id:{"type":reading_type},"Accession_run_id":{"type":reading_type}}
```

For example:

```
{"ERR3593315":{"type":"PAIRED"},"SRR7533096":{"type":"SINGLE"}}
```

The user has also to specify the name of the input file in the Snakefile (with open...).

## Running ARGfinder

In order to run ARGfinder the user should execute the following command:

```
snakemake --profile profile_argfinder
```

We remind that the execution of the pipeline requires the installation of Snakemake workflow management beforehand.

## Output

When successfully executed, ARGfinder creates a directory by the name ``` results ``` where the user can find out all the available results from all the analysis steps (results are separated in single and paired reads results). More specifically:


* ``` raw_reads ``` directory contains all the downloaded sequencing datasets + benchmarking files
* ``` trimmed_reads ``` directory contains all the trimmed sequencing datasets + benchmarking files
* ``` kma_mOTUs ``` directory contains all the alignment result files with mOTUs database.
	* res file = A result overview giving the most common statistics for each mapped template.
	* mapstat file = 
	* fsa.gz = The consensus sequences drawn from the alignments (zipped)
	* mat.gz = Base counts on each position in each template (zipped)
	* vcf.gz = 
* ``` kma_panres ``` directory contains all the alignment result files with PanRes database.
	* res file = A result overview giving the most common statistics for each mapped template.
	* mapstat file = 
	* mapstat.filtered file = 
	* bam file = 
	* fsa.gz = The consensus sequences drawn from the alignments (zipped)
	* mat.gz = Base counts on each position in each template (zipped)
	* vcf.gz = 
* ``` ARGextender ``` for extracting the genomic flanking regions around ARGs
	* fasta file = fasta file with the extracted flanking sequences
	* frag.gz file = overview file that contains information on the following: Contig_seq, Number of matching ARGs, Alignment score, Start pos., End pos., Template name, Contig name
	* frag_raw.gz = Similar file with fra.gz but with all ARGs that can align to any of the contigs
	* gfa.gz =  
* ``` Mash ``` directory contains the mash sketches for each sequecning dataset

## Tips and Tricks

* In order to save space ARGfinder removes the raw and trimmed reads from each sample. 
* Apart from the regular output files, kma_panres also provides a BAM file.
* ARGextender will run until there is nothing more to extend.
* The user should create the ``` logs ``` directory in the main directory.
* We have noticed that the enaBrowserTools conda package is a bit unstable and from time to time it is not working properly. In that case we suggest that you clone their repository in the ``` prerequisites ``` directory and use the code they provide (You can find the repository here <a href="https://github.com/enasequence/enaBrowserTools/">here</a>)

## Citation