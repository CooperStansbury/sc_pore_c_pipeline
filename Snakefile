import pandas as pd
import yaml
from pathlib import Path
from utils import snakemake_utils

# Define config
BASE_DIR = Path(workflow.basedir)
configfile: BASE_DIR/"config/config.yaml"

# get the list of basenames from the input directory
input_path = f"{BASE_DIR}/{config['reads']}"
SAMPLES = snakemake_utils.get_sample_list(input_path)


rule all:
    input:
        "outputs/tables/raw_alignment_table.csv",
        "outputs/tables/digested_fragments_table.csv",
        "outputs/tables/alignment_table.csv"
        
        
rule bwa_map:
    input:
        refgenome=f"{config['reference']}",
        reads=config['reads'] + "/{sample}.fastq"
    output:
        bam="outputs/mapped/{sample}.bam"
    log:
        "outputs/logs/{sample}.log"
    threads:
        8
    shell:
        "bwa mem -t {threads} {input.refgenome} {input.reads} "
        " | samtools view -Sb -> {output.bam} 2>{log}"
        

rule samtools_merge:
    input:
        expand("outputs/mapped/{sample}.bam", sample=SAMPLES)
    output:
        "outputs/merged.bam"
    threads:  
        8     
    wrapper:
        "0.77.0/bio/samtools/merge"


rule reformat_bam:
    input:
        "outputs/merged.bam"
    output:
        "outputs/merged_reformatted.bam"
    shell:
        "cat {input} | python3 scripts/reformat_bam.py > {output} "


rule samtools_sort:
    input:
        "outputs/merged_reformatted.bam"
    output:
        "outputs/merged_sorted.bam"
    shell:
        "samtools sort -T {input} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        "outputs/merged_sorted.bam"
    output:
        "outputs/merged_sorted.bam.bai"
    shell:
        "samtools index {input}"
        
        
rule create_table:
    input:
        "outputs/merged_sorted.bam"
    output:
        "outputs/tables/raw_alignment_table.csv"
    shell:
        "cat {input} | python3 scripts/create_table.py > {output}"
        
        
rule virtual_digest:
    input:
        f"{config['reference']}"
    output:
        "outputs/tables/digested_fragments_table.csv"
    shell:
        "python3 scripts/virtual_digest.py {input} > {output}"

        
rule assign_fragments:
    input:
        alignment_table="outputs/tables/raw_alignment_table.csv",
        fragments_table="outputs/tables/digested_fragments_table.csv"
    output:
        "outputs/tables/alignment_table.csv"
    shell:
        "python3 scripts/assign_fragments.py {input.alignment_table} {input.fragments_table} > {output}"

