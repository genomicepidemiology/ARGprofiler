import json
import glob
import os

include: "rules/analysis_paired_read.smk"
include: "rules/analysis_single_read.smk"
include: "rules/analysis_paired_read_local.smk"
include: "rules/analysis_single_read_local.smk"

#Directories of local raw reads
path_paired_raw_reads="local_reads/paired_reads/*.gz"
path_single_raw_reads="local_reads/single_reads/*.gz"

local_paired_sample_set=set()
local_single_sample=[]

for p in glob.glob(path_paired_raw_reads):
	sample_name = os.path.splitext(os.path.basename(p))[0].split('_')[0]
	local_paired_sample_set.add(sample_name)

local_single_sample = [os.path.splitext(os.path.basename(p))[0].split('.')[0] for p in glob.glob(path_single_raw_reads)]

local_paired_sample = list(local_paired_sample_set)

with open("runs_batch_test.json", 'r') as f:
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
		# expand("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_1.fastq.gz", paired_reads=paired),
		# expand("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_2.fastq.gz", paired_reads=paired),
		# expand("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_check_file_raw.txt", paired_reads=paired),
		# expand("results/raw_reads/single_end/{single_reads}/{single_reads}.fastq.gz", single_reads=single),
		# expand("results/raw_reads/single_end/{single_reads}/{single_reads}_check_file_raw.txt", single_reads=single),
		# expand("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq", paired_reads=paired),
		# expand("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq", paired_reads=paired),
		# expand("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq", paired_reads=paired),
		# expand("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_check_file_trim.txt", paired_reads=paired),
		# expand("results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq", single_reads=single),
		# expand("results/trimmed_reads/single_end/{single_reads}/{single_reads}_check_file_trim.txt", single_reads=single),
		# expand("results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.res", paired_reads=paired),
		# expand("results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.mapstat", paired_reads=paired),
		# expand("results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt", paired_reads=paired),
		# expand("results/kma_mOTUs/single_end/{single_reads}/{single_reads}.res", single_reads=single),
		# expand("results/kma_mOTUs/single_end/{single_reads}/{single_reads}.mapstat", single_reads=single),
		# expand("results/kma_mOTUs/single_end/{single_reads}/{single_reads}_check_file_kma.txt", single_reads=single),
		# expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.res", paired_reads=paired),
		# expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mat.gz", paired_reads=paired),
		# expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat", paired_reads=paired),
		# expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat.filtered", paired_reads=paired),
		# expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.bam", paired_reads=paired),
		# expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt", paired_reads=paired),
		# expand("results/kma_panres/single_end/{single_reads}/{single_reads}.res", single_reads=single),
		# expand("results/kma_panres/single_end/{single_reads}/{single_reads}.mat.gz", single_reads=single),
		# expand("results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat", single_reads=single),
		# expand("results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat.filtered", single_reads=single),
		# expand("results/kma_panres/single_end/{single_reads}/{single_reads}.bam", single_reads=single),
		# expand("results/kma_panres/single_end/{single_reads}/{single_reads}_check_file_kma.txt", single_reads=single),
		# expand("results/kma_virulence/paired_end/{paired_reads}/{paired_reads}.res", paired_reads=paired),
		# expand("results/kma_virulence/paired_end/{paired_reads}/{paired_reads}.mat.gz", paired_reads=paired),
		# expand("results/kma_virulence/paired_end/{paired_reads}/{paired_reads}.mapstat", paired_reads=paired),
		# expand("results/kma_virulence/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt", paired_reads=paired),
		# expand("results/kma_virulence/single_end/{single_reads}/{single_reads}.res", single_reads=single),
		# expand("results/kma_virulence/single_end/{single_reads}/{single_reads}.mat.gz", single_reads=single),
		# expand("results/kma_virulence/single_end/{single_reads}/{single_reads}.mapstat", single_reads=single),
		# expand("results/kma_virulence/single_end/{single_reads}/{single_reads}_check_file_kma.txt", single_reads=single),
		# expand("results/mash_sketch/paired_end/{paired_reads}/{paired_reads}.trimmed.fastq.msh", paired_reads=paired),
		# expand("results/mash_sketch/paired_end/{paired_reads}/{paired_reads}_check_file_mash.txt", paired_reads=paired),
		# expand("results/mash_sketch/single_end/{single_reads}/{single_reads}.trimmed.fastq.msh", single_reads=single),
		# expand("results/mash_sketch/single_end/{single_reads}/{single_reads}_check_file_mash.txt", single_reads=single),
		# expand("results/seed_extender/paired_end/{paired_reads}/{paired_reads}.fasta.gz", paired_reads=paired),
		# expand("results/seed_extender/paired_end/{paired_reads}/{paired_reads}.gfa.gz", paired_reads=paired),
		# expand("results/seed_extender/paired_end/{paired_reads}/{paired_reads}.frag.gz", paired_reads=paired),
		# expand("results/seed_extender/paired_end/{paired_reads}/{paired_reads}.frag_raw.gz", paired_reads=paired),
		# expand("results/seed_extender/paired_end/{paired_reads}/{paired_reads}_check_file_seed.txt", paired_reads=paired),
		# expand("results/seed_extender/single_end/{single_reads}/{single_reads}.fasta.gz", single_reads=single),
		# expand("results/seed_extender/single_end/{single_reads}/{single_reads}.gfa.gz", single_reads=single),
		# expand("results/seed_extender/single_end/{single_reads}/{single_reads}.frag.gz", single_reads=single),
		# expand("results/seed_extender/single_end/{single_reads}/{single_reads}.frag_raw.gz", single_reads=single),
		# expand("results/seed_extender/single_end/{single_reads}/{single_reads}_check_file_seed.txt", single_reads=single),
		expand("results/raw_reads/paired_end/{paired_reads}/check_clean_raw.txt", paired_reads=paired),
		expand("results/trimmed_reads/paired_end/{paired_reads}/check_clean_trim.txt", paired_reads=paired),
		expand("results/raw_reads/single_end/{single_reads}/check_clean_raw.txt", single_reads=single),
		expand("results/trimmed_reads/single_end/{single_reads}/check_clean_trim.txt", single_reads=single),
		expand("results/trimmed_reads/paired_end/{local_paired_reads}/{local_paired_reads}_1.trimmed_local.fastq", local_paired_reads=local_paired_sample),
		expand("results/trimmed_reads/paired_end/{local_paired_reads}/{local_paired_reads}_2.trimmed_local.fastq", local_paired_reads=local_paired_sample),
		expand("results/trimmed_reads/paired_end/{local_paired_reads}/{local_paired_reads}_singleton_local.trimmed.fastq", local_paired_reads=local_paired_sample),
		expand("results/trimmed_reads/paired_end/{local_paired_reads}/{local_paired_reads}_check_file_local_trim.txt", local_paired_reads=local_paired_sample),
		expand("results/kma_mOTUs/paired_end/{local_paired_reads}/{local_paired_reads}_local.res", local_paired_reads=local_paired_sample),
		expand("results/kma_mOTUs/paired_end/{local_paired_reads}/{local_paired_reads}_local.mapstat", local_paired_reads=local_paired_sample),
		expand("results/kma_mOTUs/paired_end/{local_paired_reads}/{local_paired_reads}_check_file_local_kma.txt", local_paired_reads=local_paired_sample),
		expand("results/kma_panres/paired_end/{local_paired_reads}/{local_paired_reads}_local.res", local_paired_reads=local_paired_sample),
		expand("results/kma_panres/paired_end/{local_paired_reads}/{local_paired_reads}_local.mat.gz", local_paired_reads=local_paired_sample),
		expand("results/kma_panres/paired_end/{local_paired_reads}/{local_paired_reads}_local.mapstat", local_paired_reads=local_paired_sample),
		expand("results/kma_panres/paired_end/{local_paired_reads}/{local_paired_reads}_local.bam", local_paired_reads=local_paired_sample),
		expand("results/kma_panres/paired_end/{local_paired_reads}/{local_paired_reads}_local.mapstat.filtered", local_paired_reads=local_paired_sample),
		expand("results/kma_panres/paired_end/{local_paired_reads}/{local_paired_reads}_check_file_local_kma.txt", local_paired_reads=local_paired_sample),
		expand("results/kma_virulence/paired_end/{local_paired_reads}/{local_paired_reads}_local.res", local_paired_reads=local_paired_sample),
		expand("results/kma_virulence/paired_end/{local_paired_reads}/{local_paired_reads}_local.mat.gz", local_paired_reads=local_paired_sample),
		expand("results/kma_virulence/paired_end/{local_paired_reads}/{local_paired_reads}_local.mapstat", local_paired_reads=local_paired_sample),
		expand("results/kma_virulence/paired_end/{local_paired_reads}/{local_paired_reads}_check_file_local_kma.txt", local_paired_reads=local_paired_sample),
		expand("results/mash_sketch/paired_end/{local_paired_reads}/{local_paired_reads}.trimmed_local.fastq.msh", local_paired_reads=local_paired_sample),
		expand("results/mash_sketch/paired_end/{local_paired_reads}/{local_paired_reads}_check_file_local_mash.txt", local_paired_reads=local_paired_sample),
		expand("results/seed_extender/paired_end/{local_paired_reads}/{local_paired_reads}_local.fasta.gz", local_paired_reads=local_paired_sample),
		expand("results/seed_extender/paired_end/{local_paired_reads}/{local_paired_reads}_local.gfa.gz", local_paired_reads=local_paired_sample),
		expand("results/seed_extender/paired_end/{local_paired_reads}/{local_paired_reads}_local.frag.gz", local_paired_reads=local_paired_sample),
		expand("results/seed_extender/paired_end/{local_paired_reads}/{local_paired_reads}_local.frag_raw.gz", local_paired_reads=local_paired_sample),
		expand("results/seed_extender/paired_end/{local_paired_reads}/{local_paired_reads}_check_file_local_seed.txt", local_paired_reads=local_paired_sample),
		expand("results/ppr_meta/paired_end/{local_paired_reads}/{local_paired_reads}_local.csv", local_paired_reads=local_paired_sample),
		expand("results/ppr_meta/paired_end/{local_paired_reads}/{local_paired_reads}_check_file_local_ppr.txt", local_paired_reads=local_paired_sample),
		expand("results/trimmed_reads/single_end/{local_single_reads}/{local_single_reads}.trimmed_local.fastq", local_single_reads=local_single_sample),
		expand("results/trimmed_reads/single_end/{local_single_reads}/{local_single_reads}_check_file_local_trim.txt", local_single_reads=local_single_sample),
		expand("results/kma_mOTUs/single_end/{local_single_reads}/{local_single_reads}_local.res", local_single_reads=local_single_sample),
		expand("results/kma_mOTUs/single_end/{local_single_reads}/{local_single_reads}_local.mapstat", local_single_reads=local_single_sample),
		expand("results/kma_mOTUs/single_end/{local_single_reads}/{local_single_reads}_check_file_local_kma.txt", local_single_reads=local_single_sample),
		expand("results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.res", local_single_reads=local_single_sample),
		expand("results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.mat.gz", local_single_reads=local_single_sample),
		expand("results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.mapstat", local_single_reads=local_single_sample),
		expand("results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.bam", local_single_reads=local_single_sample),
		expand("results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.mapstat.filtered", local_single_reads=local_single_sample),
		expand("results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_check_file_local_kma.txt", local_single_reads=local_single_sample),
		expand("results/kma_virulence/single_end/{local_single_reads}/{local_single_reads}_local.res", local_single_reads=local_single_sample),
		expand("results/kma_virulence/single_end/{local_single_reads}/{local_single_reads}_local.mat.gz", local_single_reads=local_single_sample),
		expand("results/kma_virulence/single_end/{local_single_reads}/{local_single_reads}_local.mapstat", local_single_reads=local_single_sample),
		expand("results/kma_virulence/single_end/{local_single_reads}/{local_single_reads}_check_file_local_kma.txt", local_single_reads=local_single_sample),
		expand("results/mash_sketch/single_end/{local_single_reads}/{local_single_reads}.trimmed_local.fastq.msh", local_single_reads=local_single_sample),
		expand("results/mash_sketch/single_end/{local_single_reads}/{local_single_reads}_check_file_local_mash.txt", local_single_reads=local_single_sample),
		expand("results/seed_extender/single_end/{local_single_reads}/{local_single_reads}_local.fasta.gz", local_single_reads=local_single_sample),
		expand("results/seed_extender/single_end/{local_single_reads}/{local_single_reads}_local.gfa.gz", local_single_reads=local_single_sample),
		expand("results/seed_extender/single_end/{local_single_reads}/{local_single_reads}_local.frag.gz", local_single_reads=local_single_sample),
		expand("results/seed_extender/single_end/{local_single_reads}/{local_single_reads}_local.frag_raw.gz", local_single_reads=local_single_sample),
		expand("results/seed_extender/single_end/{local_single_reads}/{local_single_reads}_check_file_local_seed.txt", local_single_reads=local_single_sample),
		expand("results/ppr_meta/single_end/{local_single_reads}/{local_single_reads}_local.csv", local_single_reads=local_single_sample),
		expand("results/ppr_meta/single_end/{local_single_reads}/{local_single_reads}_check_file_local_ppr.txt", local_single_reads=local_single_sample)
