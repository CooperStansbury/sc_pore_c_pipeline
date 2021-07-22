import pandas as pd
import yaml
from pathlib import Path
from utils import file_ops
from utils import GL_module_def

# Define config
BASE_DIR = Path(workflow.basedir)
configfile: BASE_DIR/"config/config.yaml"

# define modules for GL
GL_module_def.load_modules()

# build the output directory 
output_path = f"{BASE_DIR}/{config['output_dir']}"
file_ops.build_output_dir(output_path)

# get the list of samples
input_path = f"{BASE_DIR}/{config['reads_dir']}"
SAMPLES = file_ops.get_samples(input_path)


# rule all:
#     input:
#         bam="outputs/alignments/test.bam"
# 

rule index_bwa:
    input:
        refgenome=f"{config['refgenome']}"
    output:
        indexed_ref=f"{config['refgenome']}"
    shell:
        "bwa index {input.refgenome}"

# rule bwa_map:
#     input:
#         refgenome=f"{config['refgenome']}",
#         reads=expand("{input_dir}/{sample}.fastq", input_dir=config['reads_dir'], sample=SAMPLES)
#     output:
#         bam="outputs/alignments/{sample}.bam",
#     log:
#         f"{output_path}/logs/run.log"
#     shell:
#         "bwa mem {input.refgenome} {input.reads} | samtools view -Sb - > {output.bam}"        
# 




# rule bwa_map:
#     input:
#         refgenome=f"{config['refgenome']}",
#         reads=expand("{input_dir}/{sample}.fastq", input_dir=config['reads_dir'], sample=SAMPLES)
#     output:
#         bam={BASE_DIR} / {config['output_dir']},
#     log:
#         f"{output_path}/logs/run.log"
#     shell:
#         "bwa mem {input.refgenome} {input.reads} "
#         " | python3 {BASE_DIR}/scripts/reformat_bam.py "
#         " | samtools sort -O bam -o {output.bam} -);"
#         " samtools index {output.bam} "