@TODO: read_names are broken here and need to be updated
rule predigest_reads:
    input:
        f"{READS}{{sample}}.fastq"
    output:
        f"{OUTPUTS}digested/{{sample}}.fastq"
    shell:
        "python3 scripts/pre_digest.py {input} {output}"


rule filter_bookends:
    input:
        f"{OUTPUTS}tables/alignment_table.csv"
    output:
        f"{OUTPUTS}tables/filtered_alignment_table.csv"
    shell:
        f"python3 scripts/filter_bookends.py {{input}} > {{output}}"
    
        
rule build_paohviz_table:
    input:
        f"{OUTPUTS}tables/filtered_alignment_table.csv"
    output:
        f"{OUTPUTS}tables/paohviz_output.csv"
    shell:
        "python3 scripts/build_paohviz.py {input} > {output}"