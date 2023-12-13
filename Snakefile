import sys
import json
import glob
import os
import re
import subprocess

include: "rules/analysis_paired_read.smk"
include: "rules/analysis_single_read.smk"
include: "rules/fetch_db.smk"

def build_inputs(jsonFile):
	#p_id = re.compile(r'.+\/(\w+[^_1])\_?\d?\..+\.gz')
	p_id = re.compile(r'.+\/(((\w+)_\d)|(\w+))\..+\.gz')

	single, paired = set(), set()
	# First extract run_ids and run_read_types from jsonFile
	with open(jsonFile, 'r') as f:
		data = json.load(f)

	for sample_id, sample_id_read_type in data.items():
		if sample_id_read_type['type'] == 'PAIRED':
			paired.add(sample_id)
		elif sample_id_read_type['type'] == 'SINGLE':
			single.add(sample_id)

	# Directories of local raw reads
	local_raw_reads=sorted(glob.glob(os.path.join("local_reads", "*.gz")))
	i = 0
	while i < len(local_raw_reads):
		local_file = local_raw_reads[i]
		match = p_id.findall(local_file)[0]
		if len(match[-2]) > 0:
			m = match[-2]
		else:
			m = match[-1]
		read_type = None
		local_file2=""
	
		if '_1' in local_file: # paired
			# check reverse is there
			local_file2 = local_file.replace('_1', '_2')
			if local_file2 == local_raw_reads[i+1]: 	
				paired.add(m)
				read_type="paired_end"
		else:
			single.add(m)
			read_type="single_end"
		
		# move files
		if read_type in ['single_end', 'paired_end']:
			destDir = os.path.join("results", "raw_reads", read_type, m)
			checkFile = os.path.join(destDir, m + '_check_file_raw.txt')
			os.makedirs(destDir, exist_ok=True)
			p = subprocess.run(f"cp -u -t {destDir} {local_file} {local_file2}", shell=True)
			if p.returncode == 0:
				subprocess.run(f"> {checkFile}", shell=True)

		if read_type == 'paired_end':
			i += 2
		else:
			i += 1		
	return list(single), list(paired)

single, paired = build_inputs("input.json")

rule all:
	input:
		#expand("prerequisites/db_panres/panres_genes.fa"),
		expand("prerequisites/db_panres/check_file_index_db_panres.txt"),
		expand("prerequisites/db_motus/check_file_index_db_mOTUs.txt"),
		#expand("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_1.fastq.gz", paired_reads=paired),
		#expand("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_2.fastq.gz", paired_reads=paired),
		#expand("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_check_file_raw.txt", paired_reads=paired),
		#expand("results/raw_reads/single_end/{single_reads}/{single_reads}.fastq.gz", single_reads=single),
		#expand("results/raw_reads/single_end/{single_reads}/{single_reads}_check_file_raw.txt", single_reads=single),
		#expand("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq", paired_reads=paired),
		#expand("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq", paired_reads=paired),
		#expand("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq", paired_reads=paired),
		#expand("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_check_file_trim.txt", paired_reads=paired),
		#expand("results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq", single_reads=single),
		#expand("results/trimmed_reads/single_end/{single_reads}/{single_reads}_check_file_trim.txt", single_reads=single),
		#expand("results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.res", paired_reads=paired),
		#expand("results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.mapstat", paired_reads=paired),
		#expand("results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt", paired_reads=paired),
		#expand("results/kma_mOTUs/single_end/{single_reads}/{single_reads}.res", single_reads=single),
		#expand("results/kma_mOTUs/single_end/{single_reads}/{single_reads}.mapstat", single_reads=single),
		#expand("results/kma_mOTUs/single_end/{single_reads}/{single_reads}_check_file_kma.txt", single_reads=single),
		#expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.res", paired_reads=paired),
		#expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mat.gz", paired_reads=paired),
		#expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat", paired_reads=paired),
		#expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat.filtered", paired_reads=paired),
		#expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}.bam", paired_reads=paired),
		#expand("results/kma_panres/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt", paired_reads=paired),
		#expand("results/kma_panres/single_end/{single_reads}/{single_reads}.res", single_reads=single),
		#expand("results/kma_panres/single_end/{single_reads}/{single_reads}.mat.gz", single_reads=single),
		#expand("results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat", single_reads=single),
		#expand("results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat.filtered", single_reads=single),
		#expand("results/kma_panres/single_end/{single_reads}/{single_reads}.bam", single_reads=single),
		#expand("results/kma_panres/single_end/{single_reads}/{single_reads}_check_file_kma.txt", single_reads=single),
		#expand("results/mash_sketch/paired_end/{paired_reads}/{paired_reads}.trimmed.fastq.msh", paired_reads=paired),
		#expand("results/mash_sketch/paired_end/{paired_reads}/{paired_reads}_check_file_mash.txt", paired_reads=paired),
		#expand("results/mash_sketch/single_end/{single_reads}/{single_reads}.trimmed.fastq.msh", single_reads=single),
		#expand("results/mash_sketch/single_end/{single_reads}/{single_reads}_check_file_mash.txt", single_reads=single),
		#expand("results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.fasta.gz", paired_reads=paired),
		#expand("results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.gfa.gz", paired_reads=paired),
		#expand("results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.frag.gz", paired_reads=paired),
		#expand("results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.frag_raw.gz", paired_reads=paired),
		#expand("results/ARG_extender/paired_end/{paired_reads}/{paired_reads}_check_file_ARG.txt", paired_reads=paired),
		#expand("results/ARG_extender/single_end/{single_reads}/{single_reads}.fasta.gz", single_reads=single),
		#expand("results/ARG_extender/single_end/{single_reads}/{single_reads}.gfa.gz", single_reads=single),
		#expand("results/ARG_extender/single_end/{single_reads}/{single_reads}.frag.gz", single_reads=single),
		#expand("results/ARG_extender/single_end/{single_reads}/{single_reads}.frag_raw.gz", single_reads=single),
		#expand("results/ARG_extender/single_end/{single_reads}/{single_reads}_check_file_ARG.txt", single_reads=single),
		expand("results/raw_reads/paired_end/{paired_reads}/check_clean_raw.txt", paired_reads=paired),
		expand("results/trimmed_reads/paired_end/{paired_reads}/check_clean_trim.txt", paired_reads=paired),
		expand("results/raw_reads/single_end/{single_reads}/check_clean_raw.txt", single_reads=single),
		expand("results/trimmed_reads/single_end/{single_reads}/check_clean_trim.txt", single_reads=single)
