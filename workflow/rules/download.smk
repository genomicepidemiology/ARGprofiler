rule download_single_end_reads:
	"""
	Downloading metagenomic raw single end reads from ENA using enaDataGet
	"""
	output:
		"results/raw_reads/single_end/{single_reads}/{single_reads}.fastq.gz",
		check_file_raw="results/raw_reads/single_end/{single_reads}/{single_reads}_check_file_raw.txt",
		out_time="results/raw_reads/single_end/{single_reads}/{single_reads}.bench",
	envmodules:
		"tools",
		"fastq-dl/2.0.4",
	conda: "../env/environment_argprofiler.yaml"
	params:
		time=config["time_path"],
		attempts=config["max_attempts"],
		scripts_dir=config["scripts"],
		update_mysql=config["update_mysql"]
	threads: 20
	log:
		"results/raw_reads/single_end/{single_reads}/{single_reads}.log"
	shell:
		"""
		{params.time} -v --output={output.out_time} fastq-dl -a {wildcards.single_reads} --silent --cpus {threads} --max-attempts {params.attempts} -o results/raw_reads/single_end/{wildcards.single_reads} > {log}
        if [ "{params.update_mysql}" = "true" ]; then
			python {params.scripts_dir}/get_time.py -f {output.out_time} --run
		fi >> {log}

		touch {output.check_file_raw}
		"""

rule download_paired_end_reads:
	"""
	Downloading metagenomic raw paired end reads from ENA using enaDataGet
	"""
	output:
		"results/raw_reads/paired_end/{paired_reads}/{paired_reads}_1.fastq.gz",
		"results/raw_reads/paired_end/{paired_reads}/{paired_reads}_2.fastq.gz",
		check_file_raw="results/raw_reads/paired_end/{paired_reads}/{paired_reads}_check_file_raw.txt",
		out_time="results/raw_reads/paired_end/{paired_reads}/{paired_reads}.bench",
	envmodules:
		"tools",
		"fastq-dl/2.0.4",
	conda: "../env/environment_argprofiler.yaml"
	params:
		time=config["time_path"],
		attempts=config["max_attempts"],
		scripts_dir=config["scripts"],
		update_mysql=config["update_mysql"]
	threads: 20
	log:
		"results/raw_reads/paired_end/{paired_reads}/{paired_reads}.log"
	shell:
		"""
		{params.time} -v --output={output.out_time} fastq-dl -a {wildcards.paired_reads} --silent --cpus {threads} --max-attempts {params.attempts} -o results/raw_reads/paired_end/{wildcards.paired_reads} > {log}
        if [ "{params.update_mysql}" = "true" ]; then
			python {params.scripts_dir}/get_time.py -f {output.out_time} --run >> {log}
		fi
		touch {output.check_file_raw}
		"""

