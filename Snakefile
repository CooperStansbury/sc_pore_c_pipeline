import pandas as pd
import yaml
from pathlib import Path
from utils import snakemake_utils

# Define config and input and output directories
BASE_DIR = Path(workflow.basedir)
configfile: BASE_DIR/"config/config.yaml"

OUTPUTS = f"{config['output']}"
READS = f"{config['reads']}"
REFERENCE = f"{config['reference']}"
ASSEMBLY = f"{config['assembly']}"


# get the list of basenames from the input directory
SAMPLES = snakemake_utils.get_sample_list(READS)


rule all:
    input:
        f"{OUTPUTS}merged_sorted.bam",
        f"{OUTPUTS}config.yaml",
        f"{OUTPUTS}stats/merged_sorted_stats.txt",
        f"{OUTPUTS}stats/coverage_table.txt",
        f"{OUTPUTS}stats/samtools_coverage.txt",
        f"{OUTPUTS}stats/custom_coverage_table.csv",
        f"{OUTPUTS}tables/raw_alignment_table.csv",
        f"{OUTPUTS}tables/digested_fragments_table.csv",
        f"{OUTPUTS}tables/alignment_table_mapped.csv",
        f"{OUTPUTS}tables/alignment_table.csv",
        f"{OUTPUTS}tables/filtered_alignment_table.csv",
        f"{OUTPUTS}tables/paohviz_output.csv",
        f"{OUTPUTS}tables/incidence_table.csv",
        
        
rule copy_config:
    input:
        "config/config.yaml"
    output:
        f"{OUTPUTS}config.yaml"
    shell:
        "cp {input} {output}"

# @TODO: read_names are broken here and need to be updated
# rule predigest_reads:
#     input:
#         f"{READS}{{sample}}.fastq"
#     output:
#         f"{OUTPUTS}digested/{{sample}}.fastq"
#     shell:
#         "python3 scripts/pre_digest.py {input} {output}"
# 

rule aligner:
    input:
        refgenome=REFERENCE,
        reads=f"{READS}{{sample}}.fastq"
        # reads=f"{OUTPUTS}digested/{{sample}}.fastq" # for pre-digested reads 
    output:
        f"{OUTPUTS}mapped/{{sample}}.bam"
    threads:
        config['threads']
    shell:
        f"{config['aligner']} -t {{threads}} {{input.refgenome}} {{input.reads}} "
        " | samtools view -Sb -> {output}"
        

rule samtools_merge:
    input:
        expand(f"{OUTPUTS}mapped/{{sample}}.bam", sample=SAMPLES)
    output:
        f"{OUTPUTS}merged.bam"
    threads:  
        config['threads']     
    wrapper:
        "0.77.0/bio/samtools/merge"


rule reformat_bam:
    input:
        f"{OUTPUTS}merged.bam"
    output:
        f"{OUTPUTS}merged_reformatted.bam"
    shell:
        "cat {input} | python3 scripts/reformat_bam.py > {output} "


rule samtools_sort:
    input:
        f"{OUTPUTS}merged_reformatted.bam"
    output:
        f"{OUTPUTS}merged_sorted.bam"
    shell:
        "samtools sort -T {input} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        f"{OUTPUTS}merged_sorted.bam"
    output:
        f"{OUTPUTS}merged_sorted.bam.bai"
    shell:
        "samtools index {input}"
        

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
        "python3 scripts/map_assembly.py {input.align} {input.assembly} > {output}"
        
        
rule filter_bookends:
    input:
        f"{OUTPUTS}tables/alignment_table.csv"
    output:
        f"{OUTPUTS}tables/filtered_alignment_table.csv"
    shell:
        f"python3 scripts/filter_bookends.py {{input}} {config['chromosome']} > {{output}}"
    
        
rule build_paohviz_table:
    input:
        f"{OUTPUTS}tables/filtered_alignment_table.csv"
    output:
        f"{OUTPUTS}tables/paohviz_output.csv"
    shell:
        "python3 scripts/build_paohviz.py {input} > {output}"
        
        
rule build_incidence:
    input:
        f"{OUTPUTS}tables/filtered_alignment_table.csv"
    output:
        f"{OUTPUTS}tables/incidence_table.csv"
    shell:
        "python3 scripts/build_incidence_matrix.py {input} > {output}"