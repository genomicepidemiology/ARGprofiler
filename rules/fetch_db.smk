rule fetch_db:
    output:
        "prerequisites/panres/panres_genes.fa"
    params:
        zenodo_url="https://zenodo.org/records/10091602/files/panres_genes.fa"
    shell:
        """
        wget {params.zenodo_url} -P prerequisites/panres
        """

rule index_db:
    input:
        "prerequisites/panres/panres_genes.fa"
    output:
        check_file_index="prerequisites/panres/check_file_index_db.txt"
    envmodules:
        "tools",
        "kma/1.4.12a"
    shell:
        """
        cd prerequisites/panres
        sed "1d" panres_genes.fa |kma index -i -- -o pan_db
        cd ../../
        touch {output.check_file_index}
        """