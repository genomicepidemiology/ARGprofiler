rule trim_single_end_reads_local:
	"""
	Adapter trimming of raw single end local reads using fastp
	"""
	input:
		"local_reads/single_reads/{local_single_reads}.fastq.gz"
	output:
		"results/trimmed_reads/single_end/{local_single_reads}/{local_single_reads}.trimmed_local.fastq",
		check_file_trim="results/trimmed_reads/single_end/{local_single_reads}/{local_single_reads}_check_file_local_trim.txt"
	params:
		overlap_diff_limit="1",
		average_qual="20",
		length_required="50",
		cut_tail="--cut_tail",
		h="results/trimmed_reads/single_end/{local_single_reads}/{local_single_reads}.html",
		j="results/trimmed_reads/single_end/{local_single_reads}/{local_single_reads}.json"
	envmodules:
		"tools",
		"fastp/0.23.2",
		"mariadb/10.4.17",
		"mariadb-connector-c/3.3.2"
	shell:
		"""
		/usr/bin/time -v --output=results/trimmed_reads/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.bench fastp -i {input} -o {output} --overlap_diff_limit {params.overlap_diff_limit} --average_qual {params.average_qual} --length_required {params.length_required} {params.cut_tail} -h {params.h} -w 8 -j {params.j}
		touch {output.check_file_trim}
		"""

rule kma_single_end_reads_mOTUs_local:
	"""
	Mapping raw single local reads for identifying AMR using KMA with mOTUs db
	"""
	input: 
		"results/trimmed_reads/single_end/{local_single_reads}/{local_single_reads}.trimmed_local.fastq"
	output:
		"results/kma_mOTUs/single_end/{local_single_reads}/{local_single_reads}_local.res",
		"results/kma_mOTUs/single_end/{local_single_reads}/{local_single_reads}_local.mapstat",
		check_file_kma_mOTUs="results/kma_mOTUs/single_end/{local_single_reads}/{local_single_reads}_check_file_local_kma.txt"
	params:
		db="/home/databases/metagenomics/db/mOTUs_20221205/db_mOTU_20221205",
		outdir="results/kma_mOTUs/single_end/{local_single_reads}/{local_single_reads}",
		kma_params="-mem_mode -ef -1t1 -apm p -oa -matrix",
		new_out1="results/kma_mOTUs/single_end/{local_single_reads}/{local_single_reads}_local.res",
		new_out2="results/kma_mOTUs/single_end/{local_single_reads}/{local_single_reads}_local.mapstat"
	envmodules:
		"tools",
		"kma/1.4.12a",
		"mariadb/10.4.17",
		"mariadb-connector-c/3.3.2"
	shell:
		"""
		/usr/bin/time -v --output=results/kma_mOTUs/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.bench kma -i {input} -o {params.outdir} -t_db {params.db} {params.kma_params}
		rm results/kma_mOTUs/single_end/{wildcards.local_single_reads}/*.aln
		mv results/kma_mOTUs/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.res {params.new_out1}
		mv results/kma_mOTUs/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.mapstat {params.new_out2}
		gzip -f results/kma_mOTUs/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.fsa
		touch {output.check_file_kma_mOTUs}
		"""

rule kma_single_end_reads_panRes_local:
	"""
	Mapping raw single local reads for identifying AMR using KMA with panres db
	"""
	input: 
		"results/trimmed_reads/single_end/{local_single_reads}/{local_single_reads}.trimmed_local.fastq"
	output:
		"results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.res",
		"results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.mat.gz",
		"results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.mapstat",
		"results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.bam",
		"results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.mapstat.filtered",
		check_file_kma_panres="results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_check_file_local_kma.txt"
	params:
		db="/home/databases/metagenomics/db/panres_20230420/panres_20230420",
		outdir="results/kma_panres/single_end/{local_single_reads}/{local_single_reads}",
		kma_params="-ef -1t1 -nf -vcf -sam -matrix",
		mapstat="results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.mapstat",
		mapstat_filtered="results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.mapstat.filtered",
		mapstat_table="prerequisites/mapstat_filtering/pan_master_gene_tbl.tsv",
		new_out1="results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.res",
		new_out2="results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.mat.gz",
		new_out3="results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.mapstat",
		new_out4="results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.bam",
	envmodules:
		"tools",
		"kma/1.4.12a",
		"samtools/1.16",
		"gcc/9.4.0",
		"intel/perflibs/64/2020_update2",
		"R/4.3.0",
		"mariadb/10.4.17",
		"mariadb-connector-c/3.3.2"
	shell:
		"""
		/usr/bin/time -v --output=results/kma_panres/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.bench kma -i {input} -o {params.outdir} -t_db {params.db} {params.kma_params} |samtools fixmate -m - -|samtools view -u -bh -F 4|samtools sort -o results/kma_panres/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.bam
		rm results/kma_panres/single_end/{wildcards.local_single_reads}/*.aln
		mv results/kma_panres/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.res {params.new_out1}
		mv results/kma_panres/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.mat.gz {params.new_out2}
		mv results/kma_panres/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.mapstat {params.new_out3}
		mv results/kma_panres/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.bam {params.new_out4}
		gzip -f results/kma_panres/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.fsa
		Rscript prerequisites/mapstat_filtering/mapstatFilters.R -i {params.mapstat} -o {params.mapstat_filtered} -r {params.mapstat_table}
		touch {output.check_file_kma_panres}
		"""

rule kma_single_end_reads_virulence_finder_local:
	"""
	Mapping raw single local reads for identifying AMR using KMA with virulence_finder db
	"""
	input: 
		"results/trimmed_reads/single_end/{local_single_reads}/{local_single_reads}.trimmed_local.fastq"
	output:
		"results/kma_virulence/single_end/{local_single_reads}/{local_single_reads}_local.res",
		"results/kma_virulence/single_end/{local_single_reads}/{local_single_reads}_local.mat.gz",
		"results/kma_virulence/single_end/{local_single_reads}/{local_single_reads}_local.mapstat",
		check_file_kma_virulence="results/kma_virulence/single_end/{local_single_reads}/{local_single_reads}_check_file_local_kma.txt"
	params:
		db="prerequisites/virulence_finder_db/virulence_finder_db.fsa",
		outdir="results/kma_virulence/single_end/{local_single_reads}/{local_single_reads}",
		kma_params="-ef -1t1 -nf -vcf -matrix",
		new_out1="results/kma_virulence/single_end/{local_single_reads}/{local_single_reads}_local.res",
		new_out2="results/kma_virulence/single_end/{local_single_reads}/{local_single_reads}_local.mat.gz",
		new_out3="results/kma_virulence/single_end/{local_single_reads}/{local_single_reads}_local.mapstat"
	envmodules:
		"tools",
		"kma/1.4.12a",
		"mariadb/10.4.17",
		"mariadb-connector-c/3.3.2" 
	shell:
		"""
		/usr/bin/time -v --output=results/kma_virulence/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.bench kma -i {input} -o {params.outdir} -t_db {params.db} {params.kma_params}
		rm results/kma_virulence/single_end/{wildcards.local_single_reads}/*.aln
		mv results/kma_virulence/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.res {params.new_out1}
		mv results/kma_virulence/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.mat.gz {params.new_out2}
		mv results/kma_virulence/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.mapstat {params.new_out3}
		gzip -f results/kma_virulence/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.fsa
		touch {output.check_file_kma_virulence}
		"""   

rule mash_sketch_single_end_reads_local:
	"""
	Creation of mash sketches of single end local reads using mash
	"""
	input:
		"results/trimmed_reads/single_end/{local_single_reads}/{local_single_reads}.trimmed_local.fastq"
	output:
		out="results/mash_sketch/single_end/{local_single_reads}/{local_single_reads}.trimmed_local.fastq.msh",
		check_file_mash="results/mash_sketch/single_end/{local_single_reads}/{local_single_reads}_check_file_local_mash.txt"
	envmodules:
		"tools",
		"mash/2.3",
		"mariadb/10.4.17",
		"mariadb-connector-c/3.3.2"
	shell:
		"""
		/usr/bin/time -v --output=results/mash_sketch/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.bench mash sketch -k 31 -s 10000 -o {output.out} -r {input}
		touch {output.check_file_mash}
		"""

rule seed_extender_single_reads_local:
	"""
	Performing local seed extension of single local reads using perl script
	"""
	input:
		read_1="results/trimmed_reads/single_end/{local_single_reads}/{local_single_reads}.trimmed_local.fastq",
		panres_mapstat_filtered="results/kma_panres/single_end/{local_single_reads}/{local_single_reads}_local.mapstat.filtered",
	output:
		out_fasta="results/seed_extender/single_end/{local_single_reads}/{local_single_reads}_local.fasta.gz",
		out_gfa="results/seed_extender/single_end/{local_single_reads}/{local_single_reads}_local.gfa.gz",
		out_frag="results/seed_extender/single_end/{local_single_reads}/{local_single_reads}_local.frag.gz",
		out_frag_gz="results/seed_extender/single_end/{local_single_reads}/{local_single_reads}_local.frag_raw.gz",
		check_file_seed="results/seed_extender/single_end/{local_single_reads}/{local_single_reads}_check_file_local_seed.txt"
	params:
		seed="-1",
		temp_dir="results/seed_extender/single_end/{local_single_reads}/{local_single_reads}",
		db="/home/databases/metagenomics/db/panres_20230420/pan.fa"
	envmodules:
		"tools",
		"kma/1.4.12a",
		"anaconda3/2022.10",
		"spades/3.15.5",
		"fqgrep/0.0.3",
		"mariadb/10.4.17",
		"mariadb-connector-c/3.3.2"
	shell:
		"""
		if grep -q -v -m 1 "#" {input.panres_mapstat_filtered}; 
		then
			/usr/bin/time -v --output=results/seed_extender/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.bench perl prerequisites/seed_extender/targetAsm.pl {params.seed} {params.temp_dir} {params.db} {input.read_1}
			gzip -f results/seed_extender/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.fasta
			gzip -f results/seed_extender/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}.gfa
			touch {output.check_file_seed}
		else
			touch {output.out_fasta}
			touch {output.out_gfa}
			touch {output.out_frag}
			touch {output.out_frag_gz}
			touch {output.check_file_seed}
		fi
		"""

rule ppr_meta_single_reads_local:
	"""
	Identifying phages, chromosomes or plasmids with PPR-meta on local single end reads
	"""
	input:
		"results/seed_extender/single_end/{local_single_reads}/{local_single_reads}_local.fasta.gz"
	output:
		out1="results/ppr_meta/single_end/{local_single_reads}/{local_single_reads}_local.csv",
		check_file_ppr="results/ppr_meta/single_end/{local_single_reads}/{local_single_reads}_check_file_local_ppr.txt"
	params:
		threshold="0.7"
	envmodules:
		"tools",
		"mcr/R2018a",
		"ppr-meta/20210413"
	shell:
		"""
		if grep -q . {input}; then
			gzip -d {input}
			cd results/ppr_meta/single_end/{wildcards.local_single_reads}
   			ln -s ../../../../prerequisites/ppr_meta/* .
			PPR_Meta ../../../seed_extender/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}_local.fasta {wildcards.local_single_reads}_local.csv -t {params.threshold}
			rm model* predict.py
			gzip ../../../seed_extender/single_end/{wildcards.local_single_reads}/{wildcards.local_single_reads}_local.fasta 
			touch {wildcards.local_single_reads}_check_file_local_ppr.txt
    	else
    		touch {output.out1}
    		touch {output.check_file_ppr}
		fi
		"""