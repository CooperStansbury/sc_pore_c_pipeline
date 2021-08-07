rule build_incidence:
    input:
        f"{OUTPUTS}tables/alignment_table.csv"
    output:
        f"{OUTPUTS}tables/incidence_table.csv"
    shell:
        "python3 scripts/build_incidence_matrix.py {input} > {output}"