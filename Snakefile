import sys
import json
import glob
import os
import re
import subprocess
from collections import defaultdict

configfile: "config/config.yaml"

include: "rules/analysis_paired_read.smk"
include: "rules/analysis_single_read.smk"
include: "rules/fetch_db.smk"

def build_inputs(jsonFile, local_folder='local_reads'):

    # Sets to store sample IDs in
    single, paired = set(), set()

    # Open the json file
    with open(jsonFile, 'r') as jf:
        data = json.load(jf)

    # Loop through the json input and sort the file ids
    for sample_id, sample_id_read_type in data.items():
        if sample_id_read_type['type'] == 'PAIRED':
            paired.add(sample_id)
        elif sample_id_read_type['type'] == 'SINGLE':
            single.add(sample_id)
    
    # Now extract files in local_folder
    local_reads = sorted(glob.glob(os.path.join(local_folder, '*.fastq.gz')))
    single_files, paired_files = defaultdict(set), defaultdict(set)

    # Regex pattern to match the sample names and the specific suffixes
    p1 = re.compile(r'.+(R[12]).+\.fastq\.gz')
    p2 = re.compile(r'.+(_[12])\.fastq\.gz')

    # Now sort the files whether they are paired reads or not
    for fastqFile in local_reads:
        m1 = p1.match(fastqFile)
        m2 = p2.match(fastqFile)
        m = m1 if m1 else m2
        if m is not None:
            sample_id = os.path.basename(fastqFile).split(m.group(1))[0]
            sample_id = re.sub(r"([-_])$", "", sample_id)
            paired_files[sample_id].add(fastqFile)
            paired.add(sample_id)
            
        else:
            sample_id = os.path.basename(fastqFile).split(".fastq")[0]
            single_files[sample_id].add(fastqFile)
            single.add(sample_id)

    cmd_paired_template = "ln -sf {} {} && ln -sf {} {} && touch {}_check_file_raw.txt"
    #.format(os.path.realpath(fastqFile), os.path.join(dest_dir, os.path.basename(fastqFile)))
    for paired_id, paired_files in paired_files.items():
        # double check there are exactly two files stored for the id
        paired_files = list(paired_files)
        if len(paired_files) == 2:
            dest_dir = os.path.join("results", "raw_reads", "paired_end", paired_id)
            os.makedirs(dest_dir, exist_ok=True)
            cmd = cmd_paired_template.format(
                os.path.realpath(paired_files[0]),
                os.path.join(dest_dir, os.path.basename(paired_files[0])),
                os.path.realpath(paired_files[1]),
                os.path.join(dest_dir, os.path.basename(paired_files[1])),
                os.path.join(dest_dir, paired_id)
                
            )
            subprocess.run(cmd, shell=True)
        elif len(paired_files) == 1:
            paired.remove(paired_id)
            single.add(paired_id)
            single_files[paired_id].add(paired_files[0])

    # move single end files
    cmd_single_template = "ln -sf {} {} && touch {}_check_file_raw.txt"
    for single_id, single_files in single_files.items():
        single_files = list(single_files)
        if len(single_files) == 1:
            dest_dir = os.path.join("results", "raw_reads", "single_end", single_id)
            os.makedirs(dest_dir, exist_ok=True)
            cmd = cmd_single_template.format(
                os.path.realpath(single_files[0]),
                os.path.join(dest_dir, os.path.basename(single_files[0])),
                os.path.join(dest_dir, single_id)
            )
            subprocess.run(cmd, shell=True)
    
    return list(single), list(paired)

single, paired = build_inputs("input.json")

rule all:
	input:
		#expand("prerequisites/db_panres/panres_genes.fa"),
		expand("prerequisites/db_panres/check_file_index_db_panres.txt"),
		expand("prerequisites/db_panres/panres_lengths.tsv"),
		# expand("prerequisites/db_motus/check_file_index_db_mOTUs.txt"),
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
