rule download_single_end_reads:
	"""
	Downloading metagenomic raw single end reads from ENA using enaDataGet
	"""
	output:
		"results/raw_reads/single_end/{single_reads}/{single_reads}.fastq.gz",
		check_file_raw="results/raw_reads/single_end/{single_reads}/{single_reads}_check_file_raw.txt"
	envmodules:
		"tools",
		"fastq-dl/2.0.4",
	conda: "../env/environment_argprofiler.yaml"
	params:
		time=config["time_path"],
		attempts=config["max_attempts"]
	threads: 20
	log:
		"results/raw_reads/single_end/{single_reads}/{single_reads}.log"
	shell:
		"""
		{params.time} -v --output=results/raw_reads/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bench fastq-dl -a {wildcards.single_reads} --silent --cpus {threads} --max-attempts {params.attempts} -o results/raw_reads/single_end/{wildcards.single_reads} > {log}
		touch {output.check_file_raw}
		"""

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

rule kma_single_end_reads_mOTUs:
	"""
	Mapping raw single reads for identifying AMR using KMA with mOTUs db
	"""
	input: 
		read=ancient("results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq"),
		check_file_db_mOTUs="prerequisites/db_motus/check_file_index_db_mOTUs.txt"
	output:
		"results/kma_mOTUs/single_end/{single_reads}/{single_reads}.res",
		"results/kma_mOTUs/single_end/{single_reads}/{single_reads}.mapstat",
		check_file_kma_mOTUs="results/kma_mOTUs/single_end/{single_reads}/{single_reads}_check_file_kma.txt"
	params:
		db="prerequisites/db_motus/db_mOTUs",
		outdir="results/kma_mOTUs/single_end/{single_reads}/{single_reads}",
		kma_params=config["kma_params_motus"],
		time=config["time_path"]
	envmodules:
		"tools",
		"kma/1.4.12a",
	conda: "../env/environment_argprofiler.yaml"
	threads: 20
	log:
		"results/kma_mOTUs/single_end/{single_reads}/{single_reads}.log"
	shell:
		"""
		{params.time} -v --output=results/kma_mOTUs/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bench kma -i {input.read} -o {params.outdir} -t_db {params.db} {params.kma_params} -t {threads} 2>> {log}
		rm results/kma_mOTUs/single_end/{wildcards.single_reads}/*.aln
		gzip -f results/kma_mOTUs/single_end/{wildcards.single_reads}/{wildcards.single_reads}.fsa
		touch {output.check_file_kma_mOTUs}
		"""

rule kma_single_end_reads_panRes:
	"""
	Mapping raw single reads for identifying AMR using KMA with panres db
	"""
	input: 
		read=ancient("results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq"),
		check_file_db="prerequisites/db_panres/check_file_index_db_panres.txt"
	output:
		"results/kma_panres/single_end/{single_reads}/{single_reads}.res",
		"results/kma_panres/single_end/{single_reads}/{single_reads}.mat.gz",
		"results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat",
		"results/kma_panres/single_end/{single_reads}/{single_reads}.bam",
		"results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat.filtered",
		check_file_kma_panres="results/kma_panres/single_end/{single_reads}/{single_reads}_check_file_kma.txt"
	params:
		db="prerequisites/db_panres/panres",
		outdir="results/kma_panres/single_end/{single_reads}/{single_reads}",
		kma_params=config["kma_params_panres"],
		mapstat="results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat",
		mapstat_filtered="results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat.filtered",
		mapstat_table="prerequisites/db_panres/panres_lengths.tsv",
		time=config["time_path"]
	envmodules:
		"tools",
		"kma/1.4.12a",
		"samtools/1.16",
		"gcc/9.4.0",
		"intel/perflibs/64/2020_update2",
		"R/4.3.0",
	conda: "../env/environment_argprofiler.yaml"
	threads: 2
	log: 
		"results/kma_panres/single_end/{single_reads}/{single_reads}.log"
	shell:
		"""
		{params.time} -v --output=results/kma_panres/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bench kma -i {input.read} -o {params.outdir} -t_db {params.db} {params.kma_params} -t {threads} 2> {log} |samtools fixmate -m - -|samtools view -u -bh -F 4|samtools sort -o results/kma_panres/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bam 2> {log}
		rm results/kma_panres/single_end/{wildcards.single_reads}/*.aln
		gzip -f results/kma_panres/single_end/{wildcards.single_reads}/{wildcards.single_reads}.fsa
		Rscript prerequisites/mapstat_filtering/mapstatFilters.R -i {params.mapstat} -o {params.mapstat_filtered} -r {params.mapstat_table}
		touch {output.check_file_kma_panres}
		"""

rule mash_sketch_single_end_reads:
	"""
	Creation of mash sketches of single end reads using mash
	"""
	input:
		ancient("results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq")
	output:
		out="results/mash_sketch/single_end/{single_reads}/{single_reads}.trimmed.fastq.msh",
		check_file_mash="results/mash_sketch/single_end/{single_reads}/{single_reads}_check_file_mash.txt"
	envmodules:
		"tools",
		"mash/2.3",
	conda: "../env/environment_argprofiler.yaml"
	params:
		time=config["time_path"],
		k=config["mash_k"],
		s=config["mash_s"]
	threads: 20
	log:
		"results/mash_sketch/single_end/{single_reads}/{single_reads}.log"
	shell:
		"""
		{params.time} -v --output=results/mash_sketch/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bench mash sketch -k {params.k} -s {params.s} -o {output.out} -r {input} -p {threads} 2> {log}
		touch {output.check_file_mash}
		"""

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
		check_file_ARG="results/ARG_extender/single_end/{single_reads}/{single_reads}_check_file_ARG.txt"
	params:
		ARG="-1",
		temp_dir="results/ARG_extender/single_end/{single_reads}/{single_reads}",
		db="prerequisites/db_panres/panres_genes.fa",
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
		"results/ARG_extender/single_end/{single_reads}/{single_reads}.log"
	shell:
		"""
		if grep -q -v -m 1 "#" {input.panres_mapstat_filtered}; 
		then
			echo "running argextender" > {log} 
			{params.time} -v --output=results/ARG_extender/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bench perl prerequisites/ARGextender/targetAsm.pl {params.ARG} {threads} {params.temp_dir} {params.db} {input.read_1} 2>> {log}
			gzip -f results/ARG_extender/single_end/{wildcards.single_reads}/{wildcards.single_reads}.fasta 2>> {log}
			gzip -f results/ARG_extender/single_end/{wildcards.single_reads}/{wildcards.single_reads}.gfa 2>> {log}
			touch {output.check_file_ARG}
		else
			echo "not running argextender" > {log}
			touch {output.out_fasta}
			touch {output.out_gfa}
			touch {output.out_frag}
			touch {output.out_frag_gz}
			touch {output.check_file_ARG}
		fi
		"""

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
