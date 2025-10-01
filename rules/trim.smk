
rule trim_single_end_reads:
	"""
	Adapter trimming of raw paired end reads using fastp
	"""
	input:
		ancient("results/raw_reads/single_end/{single_reads}/{single_reads}.fastq.gz")
	output:
		"results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq",
		check_file_trim="results/trimmed_reads/single_end/{single_reads}/{single_reads}_check_file_trim.txt"
	params:
		overlap_diff_limit=config["overlap_diff_limit"],
		average_qual=config["average_qual"],
		length_required=config["length_required"],
		cut_tail=config["cut_tail"],
		h="results/trimmed_reads/single_end/{single_reads}/{single_reads}.html",
		j="results/trimmed_reads/single_end/{single_reads}/{single_reads}.json",
		time=config["time_path"]
	envmodules:
		"tools",
		"fastp/0.23.2",
	conda: "../env/environment_argprofiler.yaml"
	threads: 8
	log:
		"results/trimmed_reads/single_end/{single_reads}/{single_reads}.log"
	shell:
		"""
		{params.time} -v --output=results/trimmed_reads/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bench fastp -i {input} -o {output} --overlap_diff_limit {params.overlap_diff_limit} --average_qual {params.average_qual} --length_required {params.length_required} {params.cut_tail} -h {params.h} -w {threads} -j {params.j} 2> {log}
		touch {output.check_file_trim}
		"""

rule trim_paired_end_reads:
	"""
	Adapter trimming of raw paired end reads using fastp
	"""
	input:
		in1=ancient("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_1.fastq.gz"),
		in2=ancient("results/raw_reads/paired_end/{paired_reads}/{paired_reads}_2.fastq.gz")
	output:
		out1="results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq",
		out2="results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq",
		singleton="results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq",
		check_file_trim="results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_check_file_trim.txt"
	params:
		overlap_diff_limit=config["overlap_diff_limit"],
		average_qual=config["average_qual"],
		length_required=config["length_required"],
		cut_tail=config["cut_tail"],
		out_merge="results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_merged.trimmed.fastq",
		h="results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}.html",
		j="results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}.json",
		time=config["time_path"]
	envmodules:
		"tools",
		"fastp/0.23.2",
	conda: "../env/environment_argprofiler.yaml"
	threads: 8
	log:
		"results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}.log"
	shell:
		"""
		{params.time} -v --output=results/trimmed_reads/paired_end/{wildcards.paired_reads}/{wildcards.paired_reads}.bench fastp -i {input.in1} -I {input.in2} -o {output.out1} -O {output.out2} --merge --merged_out {params.out_merge} --unpaired1 {output.singleton} --unpaired2 {output.singleton} --overlap_diff_limit {params.overlap_diff_limit} --average_qual {params.average_qual} --length_required {params.length_required} {params.cut_tail} -h {params.h} -w {threads} -j {params.j} 2> {log}
		cat {params.out_merge} >> {output.singleton} 2>> {log}
		rm {params.out_merge}
		touch {output.check_file_trim}
		"""
