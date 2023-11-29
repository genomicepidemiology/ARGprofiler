rule fetch_db_panres:
    output:
        "prerequisites/panres/panres_genes.fa"
    params:
        zenodo_url="https://zenodo.org/records/10091602/files/panres_genes.fa"
    shell:
        """
        /usr/bin/time -v --output=prerequisites/panres/fetch_panres.bench wget {params.zenodo_url} -P prerequisites/panres
        """

rule index_db_panres:
    input:
        "prerequisites/panres/panres_genes.fa"
    output:
        check_file_index="prerequisites/panres/check_file_index_db_panres.txt"
    envmodules:
        "tools",
        "kma/1.4.12a"
    shell:
        """
        cd prerequisites/panres
        /usr/bin/time -v --output=index_panres.bench sed "1d" panres_genes.fa |kma index -i -- -o pan_db
        cd ../../
        touch {output.check_file_index}
        """

rule fetch_db_mOTUs:
    output:
        "prerequisites/mOTUs/db_mOTU_v3.0.1.tar.gz"
    params:
        zenodo_url="https://zenodo.org/records/5140350/files/db_mOTU_v3.0.1.tar.gz"
    shell:
        """
        /usr/bin/time -v --output=prerequisites/mOTUs/fetch_mOTUs.bench wget {params.zenodo_url} -P prerequisites/mOTUs
        """

rule index_db_mOTUs:
    input:
        "prerequisites/mOTUs/db_mOTU_v3.0.1.tar.gz"
    output:
        check_file_index="prerequisites/mOTUs/check_file_index_db_mOTUs.txt"
    envmodules:
        "tools",
        "kma/1.4.12a"
    shell:
        """
        cd prerequisites/mOTUs
        tar -xf db_mOTU_v3.0.1.tar.gz
        /usr/bin/time -v --output=index_mOTUs.bench kma index -i db_mOTU/db_mOTU_DB_CEN.fasta -o mOTUs_db
        cd ../../
        touch {output.check_file_index}
        """