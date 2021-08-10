rule aligner:
    input:
        refgenome=REFERENCE,
        reads=f"{READS}{{sample}}{config['fastq_filetype']}"
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
        