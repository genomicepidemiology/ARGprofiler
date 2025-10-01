
rule cleanup_single_end_reads:
	"""
	Removing unwanted files
	"""
	input:
		check_file_raw="results/raw_reads/single_end/{single_reads}/{single_reads}_check_file_raw.txt",
		check_file_trim="results/trimmed_reads/single_end/{single_reads}/{single_reads}_check_file_trim.txt",
		check_file_kma_mOTUs="results/kma_mOTUs/single_end/{single_reads}/{single_reads}_check_file_kma.txt",
		check_file_kma_panres="results/kma_panres/single_end/{single_reads}/{single_reads}_check_file_kma.txt",
		check_file_mash="results/mash_sketch/single_end/{single_reads}/{single_reads}_check_file_mash.txt",
		check_file_ARG="results/ARG_extender/single_end/{single_reads}/{single_reads}_check_file_ARG.txt"
	output:
		check_file_clean_final1="results/raw_reads/single_end/{single_reads}/check_clean_raw.txt",
		check_file_clean_final2="results/trimmed_reads/single_end/{single_reads}/check_clean_trim.txt"
	params:
		check_file_clean_raw1="results/raw_reads/single_end/{single_reads}/{single_reads}.fastq.gz",
		check_file_clean_trim1="results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq",
		time=config["time_path"]
	threads: 1
	log:
		"results/raw_reads/single_end/{single_reads}/clean.log"
	shell:
		"""
		cd results/raw_reads/single_end/{wildcards.single_reads}
		rm *.gz
		cd ../../../../
		touch {params.check_file_clean_raw1}
		touch {output.check_file_clean_final1}
		cd results/trimmed_reads/single_end/{wildcards.single_reads}
		rm *.fastq
		cd ../../../../
		touch {params.check_file_clean_trim1}
		touch {output.check_file_clean_final2}
		echo "Finished cleaning up." > {log}
		"""

rule cleanup_paired_end_reads:
	"""
	Removing unwanted files
	"""
	input:
		check_file_raw="results/raw_reads/paired_end/{paired_reads}/{paired_reads}_check_file_raw.txt",
		check_file_trim="results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_check_file_trim.txt",
		check_file_kma_mOTUs="results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt",
		check_file_kma_panres="results/kma_panres/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt",
		check_file_mash="results/mash_sketch/paired_end/{paired_reads}/{paired_reads}_check_file_mash.txt",
		check_file_ARG="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}_check_file_ARG.txt"
	output:
		check_file_clean_final1="results/raw_reads/paired_end/{paired_reads}/check_clean_raw.txt",
		check_file_clean_final2="results/trimmed_reads/paired_end/{paired_reads}/check_clean_trim.txt"
	params:
		check_file_clean_raw1="results/raw_reads/paired_end/{paired_reads}/{paired_reads}_1.fastq.gz",
		check_file_clean_raw2="results/raw_reads/paired_end/{paired_reads}/{paired_reads}_2.fastq.gz",
		check_file_clean_trim1="results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq",
		check_file_clean_trim2="results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq",
		check_file_clean_trim3="results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq",
		time=config["time_path"]
	log: 
	shell:
		"""
		cd results/raw_reads/paired_end/{wildcards.paired_reads}
		rm *.gz
		cd ../../../../
		touch {params.check_file_clean_raw1}
		touch {params.check_file_clean_raw2}
		touch {output.check_file_clean_final1}
		cd results/trimmed_reads/paired_end/{wildcards.paired_reads}
		rm *.fastq
		cd ../../../../
		touch {params.check_file_clean_trim1}
		touch {params.check_file_clean_trim2}
		touch {params.check_file_clean_trim3}
		touch {output.check_file_clean_final2}
		"""
