
rule kma_single_end_reads_mOTUs:
	"""
	Mapping single reads for identifying Bacteria using KMA with mOTUs db
	"""
	input: 
		read=ancient("results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq"),
		check_file_db_mOTUs="prerequisites/db_motus/check_file_index_db_mOTUs.txt"
	output:
		"results/kma_mOTUs/single_end/{single_reads}/{single_reads}.res",
		out_mapstat="results/kma_mOTUs/single_end/{single_reads}/{single_reads}.mapstat",
		check_file_kma_mOTUs="results/kma_mOTUs/single_end/{single_reads}/{single_reads}_check_file_kma.txt",
		out_time="results/kma_mOTUs/single_end/{single_reads}/{single_reads}.bench",
	params:
		db="prerequisites/db_motus/db_mOTUs",
		outdir="results/kma_mOTUs/single_end/{single_reads}/{single_reads}",
		out_fsa="results/kma_mOTUs/single_end/{single_reads}/{single_reads}.fsa",
		kma_params=config["kma_params_motus"],
		time=config["time_path"],
		scripts_dir=config["scripts"]
	envmodules:
		"tools",
		"kma/1.4.12a",
	conda: "../env/environment_argprofiler.yaml"
	threads: 20
	log:
		"results/kma_mOTUs/single_end/{single_reads}/{single_reads}.log"
	shell:
		"""
		{params.time} -v --output={output.out_time} kma -i {input.read} -o {params.outdir} -t_db {params.db} {params.kma_params} -t {threads} 2>> {log}
		rm results/kma_mOTUs/single_end/{wildcards.single_reads}/*.aln
		gzip -f {params.out_fsa} 2>> {log}
        python {params.scripts_dir}/get_time.py -f {output.out_time} --run >> {log}
        python {params.scripts_dir}/get_mapstat.py -f {output.out_mapstat} --motus --pad_taxa --run >> {log}
		touch {output.check_file_kma_mOTUs}
		"""

rule kma_single_end_reads_panRes:
	"""
	Mapping single reads for identifying AMR using KMA with panres db
	"""
	input: 
		read=ancient("results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq"),
		check_file_db="prerequisites/db_panres/check_file_index_db_panres.txt"
	output:
		"results/kma_panres/single_end/{single_reads}/{single_reads}.res",
		"results/kma_panres/single_end/{single_reads}/{single_reads}.mat.gz",
		out_mapstat="results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat",
		out_bam="results/kma_panres/single_end/{single_reads}/{single_reads}.bam",
		out_mapstat_filtered="results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat.filtered",
		check_file_kma_panres="results/kma_panres/single_end/{single_reads}/{single_reads}_check_file_kma.txt",
		out_time="results/kma_panres/single_end/{single_reads}/{single_reads}.bench"
	params:
		db="prerequisites/db_panres/panres",
		outdir="results/kma_panres/single_end/{single_reads}/{single_reads}",
		kma_params=config["kma_params_panres"],
		mapstat="results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat",
		mapstat_filtered="results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat.filtered",
		mapstat_table="prerequisites/db_panres/panres_lengths.tsv",
		out_bam="results/kma_panres/single_end/{single_reads}/{single_reads}.bam",
		out_fsa="results/kma_panres/single_end/{single_reads}/{single_reads}.fsa",
		time=config["time_path"],
		scripts_dir=config["scripts"]
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
		{params.time} -v --output={output.out_time} kma -i {input.read} -o {params.outdir} -t_db {params.db} {params.kma_params} -t {threads} 2> {log} |samtools fixmate -m - -|samtools view -u -bh -F 4|samtools sort -o {output.out_bam} 2> {log}
		rm results/kma_panres/single_end/{wildcards.single_reads}/*.aln
		gzip -f {params.out_fsa}
		Rscript {params.scripts_dir}/mapstatFilters.R -i {params.mapstat} -o {params.mapstat_filtered} -r {params.mapstat_table}
        python {params.scripts_dir}/get_time.py -f {output.out_time} --run
        python {params.scripts_dir}/get_mapstat.py -f {output.out_mapstat} --run >> {log}
		touch {output.check_file_kma_panres}
		"""

rule kma_paired_end_reads_mOTUs:
	"""
	Mapping raw paired reads for identifying AMR using KMA with mOTUs db
	"""
	input: 
		read_1=ancient("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq"),
		read_2=ancient("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq"),
		read_3=ancient("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq"),
		check_file_db_mOTUs="prerequisites/db_motus/check_file_index_db_mOTUs.txt"
	output:
		"results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.res",
		out_mapstat="results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.mapstat",
		check_file_kma_mOTUs="results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt",
		out_time="results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.bench",
	params:
		db="prerequisites/db_motus/db_mOTUs",
		outdir="results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}",
		kma_params=config["kma_params_motus"],
		out_fsa="results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.fsa",
		time=config["time_path"],
		scripts_dir=config["scripts"]
	envmodules:
		"tools",
		"kma/1.4.12a",
	conda: "../env/environment_argprofiler.yaml"
	threads: 20
	log:
		"results/kma_mOTUs/paired_end/{paired_reads}/{paired_reads}.log"
	shell:
		"""
		{params.time} -v --output={output.out_time} kma -ipe {input.read_1} {input.read_2} -i {input.read_3} -o {params.outdir} -t_db {params.db} {params.kma_params} -t {threads} 2>> {log}
		rm results/kma_mOTUs/paired_end/{wildcards.paired_reads}/*.aln
		gzip -f {params.out_fsa} 2>> {log}
        python {params.scripts_dir}/get_time.py -f {output.out_time} --run
        python {params.scripts_dir}/get_mapstat.py -f {output.out_mapstat} --motus --pad_taxa --run >> {log}
		touch {output.check_file_kma_mOTUs}
		"""

rule kma_paired_end_reads_panRes:
	"""
	Mapping raw paired reads for identifying AMR using KMA with panres db
	"""
	input: 
		read_1=ancient("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_1.trimmed.fastq"),
		read_2=ancient("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_2.trimmed.fastq"),
		read_3=ancient("results/trimmed_reads/paired_end/{paired_reads}/{paired_reads}_singleton.trimmed.fastq"),
		check_file_db_panres="prerequisites/db_panres/check_file_index_db_panres.txt"
	output:
		"results/kma_panres/paired_end/{paired_reads}/{paired_reads}.res",
		"results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mat.gz",
		out_mapstat="results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat",
		out_mapstat_filtered="results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat.filtered",
		out_bam="results/kma_panres/paired_end/{paired_reads}/{paired_reads}.bam",
		check_file_kma_panres="results/kma_panres/paired_end/{paired_reads}/{paired_reads}_check_file_kma.txt",
		out_time="results/kma_panres/paired_end/{paired_reads}/{paired_reads}.bench"
	params:
		db="prerequisites/db_panres/panres",
		outdir="results/kma_panres/paired_end/{paired_reads}/{paired_reads}",
		kma_params=config["kma_params_panres"],
		mapstat="results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat",
		mapstat_filtered="results/kma_panres/paired_end/{paired_reads}/{paired_reads}.mapstat.filtered",
		out_fsa="results/kma_panres/paired_end/{paired_reads}/{paired_reads}.fsa",
		mapstat_table="prerequisites/db_panres/panres_lengths.tsv",
		time=config["time_path"],
		scripts_dir=config["scripts"]
	envmodules:
		"tools",
		"kma/1.4.12a",
		"samtools/1.16",
		"gcc/9.4.0",
		"intel/perflibs/64/2020_update2",
		"R/4.3.0",
	threads: 2
	conda: "../env/environment_argprofiler.yaml"
	log:
		"results/kma_panres/paired_end/{paired_reads}/{paired_reads}.log"
	shell:
		"""
		{params.time} -v --output={output.out_time} kma -ipe {input.read_1} {input.read_2} -i {input.read_3} -o {params.outdir} -t_db {params.db} {params.kma_params} -t {threads} 2> {log} |samtools fixmate -m - -|samtools view -u -bh -F 4|samtools sort -o {output.out_bam} 2>> {log}
		rm results/kma_panres/paired_end/{wildcards.paired_reads}/*.aln 
		gzip -f {params.out_fsa} 2>> {log}
		Rscript prerequisites/mapstat_filtering/mapstatFilters.R -i {params.mapstat} -o {params.mapstat_filtered} -r {params.mapstat_table} 2>> {log}
        python {params.scripts_dir}/get_time.py -f {output.out_time} --run
        python {params.scripts_dir}/get_mapstat.py -f {output.out_mapstat} --run >> {log}
		touch {output.check_file_kma_panres}
		"""
