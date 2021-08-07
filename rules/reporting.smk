rule get_stats:
    input:
        f"{OUTPUTS}merged_sorted.bam"
    output:
        f"{OUTPUTS}stats/merged_sorted_stats.txt"
    shell:
        "samtools stats {input} > {output}"
     

rule get_coverage_table:
    input:
        f"{OUTPUTS}merged_sorted.bam"
    output:
        f"{OUTPUTS}stats/coverage_table.txt"
    shell:
        "genomeCoverageBed -ibam {input} -bga > {output}"
 
 
rule per_xsome_coverage:
    input:
        f"{OUTPUTS}stats/coverage_table.txt"
    output:
        f"{OUTPUTS}stats/custom_coverage_table.csv"
    shell:
         "python3 scripts/per_xsome_coverage.py {input} > {output}"


rule samtools_coverage:
    input:
        f"{OUTPUTS}merged_sorted.bam"
    output:
        f"{OUTPUTS}stats/samtools_coverage.txt"
    shell:
        "samtools coverage {input} > {output}"