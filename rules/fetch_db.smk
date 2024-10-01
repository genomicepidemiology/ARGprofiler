rule fetch_db_panres_fa:
    output:
        "prerequisites/db_panres/panres_genes.fa"
    params:
        zenodo_url="https://zenodo.org/records/10091602/files/panres_genes.fa"
    threads: 1
    log:
        "prerequisites/db_panres/panres_genes.log"
    shell:
        """
        /usr/bin/time -v --output=prerequisites/db_panres/fetch_panres_fa.bench wget {params.zenodo_url} -P prerequisites/db_panres >> {log}
        sed -i '/^#/d' {output}
        """

rule fetch_db_panres_meta:
    output:
        glengths = "prerequisites/db_panres/panres_lengths.tsv"
    params:
        zenodo_url="https://zenodo.org/records/10091602/files/panres_annotations.tsv",
        meta = "prerequisites/db_panres/panres_annotations.tsv"
    threads: 1
    log:
        "prerequisites/db_panres/panres_meta.log"
    shell:
        """
        /usr/bin/time -v --output=prerequisites/db_panres/fetch_panres_meta.bench wget {params.zenodo_url} -P prerequisites/db_panres > {log}
        grep 'gene_length' {params.meta} | cut -f1,3 >> {output.glengths}
        rm {params.meta}
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
    log:
        "prerequisites/db_panres/panres_index.log"
    shell:
        """
        /usr/bin/time -v --output=index_panres.bench kma index -i {input} -o prerequisites/db_panres/panres 2> {log}
        touch {output.check_file_index}
        """

rule fetch_db_mOTUs:
    output:
        "prerequisites/db_motus/db_mOTU_v3.0.1.tar.gz"
    params:
        zenodo_url="https://zenodo.org/records/5140350/files/db_mOTU_v3.0.1.tar.gz"
    threads: 1
    log:
        "prerequisites/db_motus/db_mOTU_v3.0.1.log"
    shell:
        """
        /usr/bin/time -v --output=prerequisites/db_motus/fetch_mOTUs.bench wget {params.zenodo_url} -P prerequisites/db_motus >> {log}
        """

rule index_db_mOTUs:
    input:
        "prerequisites/db_motus/db_mOTU_v3.0.1.tar.gz"
    output:
        check_file_index="prerequisites/db_motus/check_file_index_db_mOTUs.txt"
    envmodules:
        "tools",
        "kma/1.4.12a"
    threads: 20
    log:
        "prerequisites/db_motus/index_db_mOTUs.log"
    shell:
        """
        tar -xf {input} -C prerequisites/db_motus/ db_mOTU/db_mOTU_DB_CEN.fasta > {log}
        /usr/bin/time -v --output=prerequisites/db_motus/index_mOTUs.bench kma index -i prerequisites/db_motus/db_mOTU/db_mOTU_DB_CEN.fasta -o prerequisites/db_motus/db_mOTUs 2>> {log}
        touch {output.check_file_index} 
        """
