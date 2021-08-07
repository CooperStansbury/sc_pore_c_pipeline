# @TODO: switch out for stable pore-c tool?
# conda install -c bioconda pore-c
rule create_table:
    input:
        f"{OUTPUTS}merged_sorted.bam"
    output:
        f"{OUTPUTS}tables/raw_alignment_table.csv"
    shell:
        "cat {input} | python3 scripts/create_table.py > {output}"
        
        
rule virtual_digest:
    input:
        REFERENCE
    output:
        f"{OUTPUTS}tables/digested_fragments_table.csv"
    shell:
        "python3 scripts/virtual_digest.py {input} > {output}"

        
rule assign_fragments:
    input:
        alignment_table=f"{OUTPUTS}tables/raw_alignment_table.csv",
        fragments_table=f"{OUTPUTS}tables/digested_fragments_table.csv"
    output:
        f"{OUTPUTS}tables/alignment_table_mapped.csv"
    shell:
        "python3 scripts/assign_fragments.py {input.alignment_table} {input.fragments_table} > {output}"
        

rule map_assembly:
    input:
        align=f"{OUTPUTS}tables/alignment_table_mapped.csv",
        assembly=ASSEMBLY
    output:
        f"{OUTPUTS}tables/alignment_table.csv"
    shell:
        f"python3 scripts/map_assembly.py {{input.align}} {{input.assembly}} {config['assembly_field']} > {{output}}"