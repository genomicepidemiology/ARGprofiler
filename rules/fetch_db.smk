rule fetch_db_panres:
    output:
        "prerequisites/db_panres/panres_genes.fa"
    params:
        zenodo_url="https://zenodo.org/records/10091602/files/panres_genes.fa"
    threads: 1
    shell:
        """
        /usr/bin/time -v --output=prerequisites/db_panres/fetch_panres.bench wget {params.zenodo_url} -P prerequisites/db_panres
        """

rule index_db_panres:
    input:
        "prerequisites/db_panres/panres_genes.fa"
    output:
        check_file_index="prerequisites/db_panres/check_file_index_db_panres.txt"
    envmodules:
        "tools",
        "kma/1.4.12a"
    threads: 1
    shell:
        """
        cd prerequisites/db_panres
        /usr/bin/time -v --output=index_panres.bench sed "1d" panres_genes.fa |kma index -i -- -o panres
        cd ../../
        touch {output.check_file_index}
        """

rule fetch_db_mOTUs:
    output:
        "prerequisites/db_motus/db_mOTU_v3.0.1.tar.gz"
    params:
        zenodo_url="https://zenodo.org/records/5140350/files/db_mOTU_v3.0.1.tar.gz"
    threads: 1
    shell:
        """
        /usr/bin/time -v --output=prerequisites/db_motus/fetch_mOTUs.bench wget {params.zenodo_url} -P prerequisites/db_motus
        """

rule index_db_mOTUs:
    input:
        "prerequisites/db_motus/db_mOTU_v3.0.1.tar.gz"
    output:
        check_file_index="prerequisites/db_motus/check_file_index_db_mOTUs.txt"
    envmodules:
        "tools",
        "kma/1.4.12a"
    threads: 1
    shell:
        """
        cd prerequisites/db_motus
        tar -xf db_mOTU_v3.0.1.tar.gz
        /usr/bin/time -v --output=index_mOTUs.bench kma index -i db_mOTU/db_mOTU_DB_CEN.fasta -o db_mOTUs
        cd ../../
        touch {output.check_file_index}
        """
