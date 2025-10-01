
rule ARG_extender_single_reads:
	"""
	Performing local ARG extension of paired reads using perl script
	"""
	input:
		read_1=ancient("results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq"),
		panres_mapstat_filtered="results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat.filtered"
	output:
		out_fasta="results/ARG_extender/single_end/{single_reads}/{single_reads}.fasta.gz",
		out_gfa="results/ARG_extender/single_end/{single_reads}/{single_reads}.gfa.gz",
		out_frag="results/ARG_extender/single_end/{single_reads}/{single_reads}.frag.gz",
		out_frag_gz="results/ARG_extender/single_end/{single_reads}/{single_reads}.frag_raw.gz",
		check_file_ARG="results/ARG_extender/single_end/{single_reads}/{single_reads}_check_file_ARG.txt",
	params:
		ARG="-1",
		temp_dir="results/ARG_extender/single_end/{single_reads}/{single_reads}",
		db="prerequisites/db_panres/panres_genes.fa",
		out_fasta="results/ARG_extender/single_end/{single_reads}/{single_reads}.fasta",
		out_gfa="results/ARG_extender/single_end/{single_reads}/{single_reads}.gfa",
		out_time="results/ARG_extender/single_end/{single_reads}/{single_reads}.bench",
		time=config["time_path"],
	envmodules:
		"tools",
		"kma/1.4.12a",
		"anaconda3/2022.10",
		"spades/3.15.5",
		"fqgrep/0.0.3",
	conda: "../env/environment_argprofiler.yaml"
	threads: 20
	log:
		"results/ARG_extender/single_end/{single_reads}/{single_reads}.log"
	shell:
		"""
		if grep -q -v -m 1 "#" {input.panres_mapstat_filtered}; 
		then
			echo "running argextender" > {log} 
			{params.time} -v --output={params.out_time} perl prerequisites/ARGextender/targetAsm.pl {params.ARG} {threads} {params.temp_dir} {params.db} {input.read_1} 2>> {log}
			gzip -f {params.out_fasta} 2>> {log}
			gzip -f {params.out_gfa} 2>> {log}
            python prerequisites/scripts_sql/get_time.py -f {params.out_time} --run
			touch {output.check_file_ARG}
		else
			echo "not running argextender" > {log}
			touch {output.out_fasta}
			touch {output.out_gfa}
			touch {output.out_frag}
			touch {output.out_frag_gz}
            python prerequisites/scripts_sql/job_status.py --run_accession {wildcards.single_reads} --rule ARG_extender_single_reads
			touch {output.check_file_ARG}
		fi
		"""

rule ARG_extender_paired_reads:
	"""
	Performing local ARG extension of paired reads using perl script
	"""
	input:
		read_1=ancient("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq"),
		read_2=ancient("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq"),
		read_3=ancient("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq"),
		panres_mapstat_filtered="results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat.filtered"
	output:
		out_fasta="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.fasta.gz",
		out_gfa="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.gfa.gz",
		out_frag="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.frag.gz",
		out_frag_gz="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.frag_raw.gz",
		check_file_ARG="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}_check_file_ARG.txt"
	params:
		ARG="-1",
		temp_dir="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}",
		out_fasta="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.fasta",
		out_gfa="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.gfa",
		db="prerequisites/db_panres/panres_genes.fa",
		out_time="results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.bench",
		time=config["time_path"]
	envmodules:
		"tools",
		"kma/1.4.12a",
		"anaconda3/2022.10",
		"spades/3.15.5",
		"fqgrep/0.0.3",
	conda: "../env/environment_argprofiler.yaml"
	threads: 20
	log:
		"results/ARG_extender/paired_end/{paired_reads}/{paired_reads}.log"
	shell:
		"""
		if grep -q -v -m 1 "#" {input.panres_mapstat_filtered}; 
		then
			echo "running argextender" > {log} 
			{params.time} -v --output={params.out_time} perl prerequisites/ARGextender/targetAsm.pl {params.ARG} {threads} {params.temp_dir} {params.db} {input.read_1} {input.read_2} {input.read_3} 2>> {log}
			gzip -f {params.out_fasta} 2>> {log}
			gzip -f {params.out_gfa} 2>> {log}
            python prerequisites/scripts_sql/get_time.py -f {params.out_time} --run
			touch {output.check_file_ARG}
		else
			echo "not running argextender" > {log}
			touch {output.out_fasta}
			touch {output.out_gfa}
			touch {output.out_frag}
			touch {output.out_frag_gz}
            python prerequisites/scripts_sql/job_status.py --run_accession {wildcards.paired_reads} --rule ARG_extender_paired_reads
			touch {output.check_file_ARG}
		fi
		"""
