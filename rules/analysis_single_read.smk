rule download_single_end_reads:
	"""
	Downloading metagenomic raw single end reads from ENA using enaDataGet
	"""
	output:
		"results/raw_reads/single_end/{single_reads}/{single_reads}.fastq.gz",
		check_file_raw="results/raw_reads/single_end/{single_reads}/{single_reads}_check_file_raw.txt"
	params:
		outdir="results/raw_reads/single_end",
		type="fastq"
	shell:
		"""
		/usr/bin/time -v --output=results/raw_reads/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bench python3 prerequisites/enaBrowserTools/python3/enaDataGet.py -f {params.type} -d {params.outdir} {wildcards.single_reads}
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
		overlap_diff_limit="1",
		average_qual="20",
		length_required="50",
		cut_tail="--cut_tail",
		h="results/trimmed_reads/single_end/{single_reads}/{single_reads}.html",
		j="results/trimmed_reads/single_end/{single_reads}/{single_reads}.json"
	conda:"environment_argfinder.yaml"
	shell:
		"""
		/usr/bin/time -v --output=results/trimmed_reads/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bench fastp -i {input} -o {output} --overlap_diff_limit {params.overlap_diff_limit} --average_qual {params.average_qual} --length_required {params.length_required} {params.cut_tail} -h {params.h} -j {params.j}
		touch {output.check_file_trim}
		"""

rule kma_single_end_reads_mOTUs:
	"""
	Mapping raw single reads for identifying AMR using KMA with mOTUs db
	"""
	input: 
		ancient("results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq")
	output:
		"results/kma_mOTUs/single_end/{single_reads}/{single_reads}.res",
		"results/kma_mOTUs/single_end/{single_reads}/{single_reads}.mapstat",
		check_file_kma_mOTUs="results/kma_mOTUs/single_end/{single_reads}/{single_reads}_check_file_kma.txt"
	params:
		db="prerequisites/db_motus/db_mOTUs",
		outdir="results/kma_mOTUs/single_end/{single_reads}/{single_reads}",
		kma_params="-mem_mode -ef -1t1 -apm p -oa -matrix"
	conda:"environment_argfinder.yaml"
	shell:
		"""
		/usr/bin/time -v --output=results/kma_mOTUs/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bench kma -i {input} -o {params.outdir} -t_db {params.db} {params.kma_params}
		rm results/kma_mOTUs/single_end/{wildcards.single_reads}/*.aln
		gzip -f results/kma_mOTUs/single_end/{wildcards.single_reads}/{wildcards.single_reads}.fsa
		touch {output.check_file_kma_mOTUs}
		"""

rule kma_single_end_reads_panRes:
	"""
	Mapping raw single reads for identifying AMR using KMA with panres db
	"""
	input: 
		ancient("results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq")
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
		kma_params="-ef -1t1 -nf -vcf -sam -matrix",
		mapstat="results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat",
		mapstat_filtered="results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat.filtered",
		mapstat_table="prerequisites/mapstat_filtering/pan_master_gene_tbl.tsv"
	conda:"environment_argfinder.yaml"
	shell:
		"""
		/usr/bin/time -v --output=results/kma_panres/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bench kma -i {input} -o {params.outdir} -t_db {params.db} {params.kma_params} |samtools fixmate -m - -|samtools view -u -bh -F 4|samtools sort -o results/kma_panres/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bam
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
	conda:"environment_argfinder.yaml"
	shell:
		"""
		/usr/bin/time -v --output=results/mash_sketch/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bench mash sketch -k 31 -s 10000 -o {output.out} -r {input}
		touch {output.check_file_mash}
		"""

rule arg_extender_single_reads:
	"""
	Performing local seed extension of paired reads using perl script
	"""
	input:
		read_1=ancient("results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq"),
		panres_mapstat_filtered="results/kma_panres/single_end/{single_reads}/{single_reads}.mapstat.filtered"
	output:
		out_fasta="results/argextender/single_end/{single_reads}/{single_reads}.fasta.gz",
		out_gfa="results/argextender/single_end/{single_reads}/{single_reads}.gfa.gz",
		out_frag="results/argextender/single_end/{single_reads}/{single_reads}.frag.gz",
		out_frag_gz="results/argextender/single_end/{single_reads}/{single_reads}.frag_raw.gz",
		check_file_seed="results/argextender/single_end/{single_reads}/{single_reads}_check_file_seed.txt"
	params:
		seed="-1",
		temp_dir="results/argextender/single_end/{single_reads}/{single_reads}",
		db="prerequisites/db_panres/pan.fa"
	conda:"environment_argfinder.yaml"
	shell:
		"""
		if grep -q -v -m 1 "#" {input.panres_mapstat_filtered}; 
		then
			/usr/bin/time -v --output=results/argextender/single_end/{wildcards.single_reads}/{wildcards.single_reads}.bench perl prerequisites/ARGextender/targetAsm.pl {params.seed} {params.temp_dir} {params.db} {input.read_1}
			gzip -f results/argextender/single_end/{wildcards.single_reads}/{wildcards.single_reads}.fasta
			gzip -f results/argextender/single_end/{wildcards.single_reads}/{wildcards.single_reads}.gfa
			touch {output.check_file_seed}
		else
			touch {output.out_fasta}
			touch {output.out_gfa}
			touch {output.out_frag}
			touch {output.out_frag_gz}
			touch {output.check_file_seed}
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
		check_file_seed="results/argextender/single_end/{single_reads}/{single_reads}_check_file_seed.txt"
	output:
		check_file_clean_final1="results/raw_reads/single_end/{single_reads}/check_clean_raw.txt",
		check_file_clean_final2="results/trimmed_reads/single_end/{single_reads}/check_clean_trim.txt"
	params:
		check_file_clean_raw1="results/raw_reads/single_end/{single_reads}/{single_reads}.fastq.gz",
		check_file_clean_trim1="results/trimmed_reads/single_end/{single_reads}/{single_reads}.trimmed.fastq"
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
		"""
