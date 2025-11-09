
rule mash_sketch_single_end_reads:
	"""
	Creation of mash sketches of single end reads using mash
	"""
	input:
		ancient("results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq")
	output:
		out="results/mash_sketch/single_end/{single_reads}/{single_reads}.trimmed.fastq.msh",
		check_file_mash="results/mash_sketch/single_end/{single_reads}/{single_reads}_check_file_mash.txt",
		out_time="results/mash_sketch/single_end/{single_reads}/{single_reads}.bench",
	envmodules:
		"tools",
		"mash/2.3",
	conda: "../env/environment_argprofiler.yaml"
	params:
		time=config["time_path"],
		k=config["mash_k"],
		s=config["mash_s"],
		scripts_dir=config["scripts"],
		update_mysql=config["update_mysql"]
	threads: 20
	log:
		"results/mash_sketch/single_end/{single_reads}/{single_reads}.log"
	shell:
		"""
		{params.time} -v --output={output.out_time} mash sketch -k {params.k} -s {params.s} -o {output.out} -r {input} -p {threads} 2> {log}
		if [ "{params.update_mysql}" = "true" ]; then
			python {params.scripts_dir}/get_time.py -f {output.out_time} --run
		fi
		touch {output.check_file_mash}
		"""

rule mash_sketch_paired_end_reads:
	"""
	Creation of mash sketches of paired end reads using mash
	"""
	input:
		read_1=ancient("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq"),
		read_2=ancient("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq"),
		read_3=ancient("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq")
	output:
		out="results/mash_sketch/paired_end/{paired_reads}/{paired_reads}.trimmed.fastq.msh",
		check_file_mash="results/mash_sketch/paired_end/{paired_reads}/{paired_reads}_check_file_mash.txt",
		out_time="results/mash_sketch/paired_end/{paired_reads}/{paired_reads}.bench"
	envmodules:
		"tools",
		"mash/2.3",
	conda: "../env/environment_argprofiler.yaml"
	params:
		time=config["time_path"],
		k=config["mash_k"],
		s=config["mash_s"],
		scripts_dir=config["scripts"],
		update_mysql=config["update_mysql"]
	threads: 20
	log:
		"results/mash_sketch/paired_end/{paired_reads}/{paired_reads}.log"
	shell:
		"""
		{params.time} -v --output={output.out_time} cat {input.read_1} {input.read_2} {input.read_3} | mash sketch -k {params.k} -s {params.s} -I {wildcards.paired_reads} -C Paired -r -o {output.out} -p {threads} - 2>> {log}
        if [ "{params.update_mysql}" = "true" ]; then
			python {params.scripts_dir}/get_time.py -f {output.out_time} --run
		fi
		touch {output.check_file_mash}
		"""
