rule build_incidence:
    input:
        f"{OUTPUTS}tables/alignment_table.csv"
    output:
        f"{OUTPUTS}tables/incidence_table.csv"
    shell:
        "python3 scripts/build_incidence_matrix.py {input} > {output}"


# rule build_chromesizes:
#     input:
#         refgenome=REFERENCE
#     output:
#         f"{OUTPUTS}hic_files/chrom_sizes"
#     shell:
#         "bioawk -c fastx '{print $name"\t"length($seq)}' {{input}} > {{output}}"
# 
# 
# rule build_pairs:
#     input:
#         alignment=f"{OUTPUTS}tables/alignment_table.csv",
#         assembly=ASSEMBLY
#     output:
#         f"{OUTPUTS}tables/pairs_file.pairs"
#     shell:
#         f"python3 scripts/build_pairs.py {{input.alignment}} {{input.assembly}} {config['assembly_field']} > {{output}}"