import json
import os

configfile: "rules/environment_argfinder.yaml"

include: "rules/analysis_paired_read.smk"
include: "rules/analysis_single_read.smk"

with open("run_batch.json", 'r') as f:
  data = json.load(f)

accession=[]
accession_type=[]

for sample_id in data:
	accession.append(sample_id)
for sample_id_type in data.values():
	accession_type.append(sample_id_type["type"])

merge=zip(accession,accession_type)
new_merge=dict(merge)

single=[]
paired=[]

for k,v in new_merge.items():
	if ("SINGLE") in v:
		single.append(k)
	elif ("PAIRED") in v:
		paired.append(k)

rule all:
	input:
		expand("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_1.fastq.gz", paired_reads=paired),
		expand("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_2.fastq.gz", paired_reads=paired),
		expand("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_check_file_raw.txt", paired_reads=paired),
		expand("results/raw_reads/single_end/{single_reads}/{single_reads}.fastq.gz", single_reads=single),
		expand("results/raw_reads/single_end/{single_reads}/{single_reads}_check_file_raw.txt", single_reads=single),
		expand("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq", paired_reads=paired),
		expand("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq", paired_reads=paired),
		expand("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq", paired_reads=paired),
		expand("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_check_file_trim.txt", paired_reads=paired),
		expand("results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq", single_reads=single),
		expand("results/trimmed_reads/single_end/{single_reads}/{single_reads}_check_file_trim.txt", single_reads=single),
		expand("results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.res", paired_reads=paired),
		expand("results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.mapstat", paired_reads=paired),
		expand("results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt", paired_reads=paired),
		expand("results/kma_mOTUs/single_end/{single_reads}/{single_reads}.res", single_reads=single),
		expand("results/kma_mOTUs/single_end/{single_reads}/{single_reads}.mapstat", single_reads=single),
		expand("results/kma_mOTUs/single_end/{single_reads}/{single_reads}_check_file_kma.txt", single_reads=single),
		expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.res", paired_reads=paired),
		expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mat.gz", paired_reads=paired),
		expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat", paired_reads=paired),
		expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat.filtered", paired_reads=paired),
		expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.bam", paired_reads=paired),
		expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt", paired_reads=paired),
		expand("results/kma_panres/single_end/{single_reads}/{single_reads}.res", single_reads=single),
		expand("results/kma_panres/single_end/{single_reads}/{single_reads}.mat.gz", single_reads=single),
		expand("results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat", single_reads=single),
		expand("results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat.filtered", single_reads=single),
		expand("results/kma_panres/single_end/{single_reads}/{single_reads}.bam", single_reads=single),
		expand("results/kma_panres/single_end/{single_reads}/{single_reads}_check_file_kma.txt", single_reads=single),
		expand("results/mash_sketch/paired_end/{paired_reads}/{paired_reads}.trimmed.fastq.msh", paired_reads=paired),
		expand("results/mash_sketch/paired_end/{paired_reads}/{paired_reads}_check_file_mash.txt", paired_reads=paired),
		expand("results/mash_sketch/single_end/{single_reads}/{single_reads}.trimmed.fastq.msh", single_reads=single),
		expand("results/mash_sketch/single_end/{single_reads}/{single_reads}_check_file_mash.txt", single_reads=single),
		expand("results/argextender/paired_end/{paired_reads}/{paired_reads}.fasta.gz", paired_reads=paired),
		expand("results/argextender/paired_end/{paired_reads}/{paired_reads}.gfa.gz", paired_reads=paired),
		expand("results/argextender/paired_end/{paired_reads}/{paired_reads}.frag.gz", paired_reads=paired),
		expand("results/argextender/paired_end/{paired_reads}/{paired_reads}.frag_raw.gz", paired_reads=paired),
		expand("results/argextender/paired_end/{paired_reads}/{paired_reads}_check_file_seed.txt", paired_reads=paired),
		expand("results/argextender/single_end/{single_reads}/{single_reads}.fasta.gz", single_reads=single),
		expand("results/argextender/single_end/{single_reads}/{single_reads}.gfa.gz", single_reads=single),
		expand("results/argextender/single_end/{single_reads}/{single_reads}.frag.gz", single_reads=single),
		expand("results/argextender/single_end/{single_reads}/{single_reads}.frag_raw.gz", single_reads=single),
		expand("results/argextender/single_end/{single_reads}/{single_reads}_check_file_seed.txt", single_reads=single),
		expand("results/raw_reads/paired_end/{paired_reads}/check_clean_raw.txt", paired_reads=paired),
		expand("results/trimmed_reads/paired_end/{paired_reads}/check_clean_trim.txt", paired_reads=paired),
		expand("results/raw_reads/single_end/{single_reads}/check_clean_raw.txt", single_reads=single),
		expand("results/trimmed_reads/single_end/{single_reads}/check_clean_trim.txt", single_reads=single)
