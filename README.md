# Single Cell Pore-C Pipeline

## Overview 

This tool is intended to create analysis files for single-cell pore-C reads.

## Dependencies

To run the whole pipeline make sure that the following modules are loaded in the GreatLakes environment:

- Bioinformatics
- bwa
- samtools/1.13
- bedtools2

Note that you can load these modules by running: `./build_gl_env.sh`.


## Usage

This is a [snakemake tool](https://snakemake.readthedocs.io/en/stable/index.html). It is designed to take reads after basecalling through to several analysis files. 
### Running the Pipeline

The whole pipeline may be run from inside the `sc_pore_c_pipeline` directory like so:

```
snakemake -j <n_cores>
```

Individual steps may be run like so, as long as they do not contain wildcards (alignment tools).

```
snakemake -j <n_cores> <rule_name>
```

### Inputs

Each input can be specified separately in the `config`. There are 4 input requirements:

- An assembly file. For example: [Mouse Genome Assembly GRCm39](https://www.ncbi.nlm.nih.gov/grc/mouse/data). Note that this file has been trimmed of extra commas and leading whitespace must be handled. Column names are expected: `Chromosome,Total length (bp),GenBank accession,RefSeq accession`
- An indexed reference genome. Note that this should be indexed by the alignment tool (e.g., `bwa`, `minimap2`).
- A directory of reads (`.fastq`) after basecalling
- A config file, from the `sc_pore_c_pipeline/config` location specifying the inputs and outputs.

```
inputs/
├── assembly
│   └── GRCm39_Assembly.csv
├── reads
│   ├── FAP70729_pass_barcode11_003e2571_0.fastq
│   ├── FAP70729_pass_barcode11_003e2571_10.fastq
│   ├── FAP70729_pass_barcode11_003e2571_11.fastq
│   ├── FAP70729_pass_barcode11_003e2571_12.fastq
│   ├── FAP70729_pass_barcode11_003e2571_13.fastq
│   ├── FAP70729_pass_barcode11_003e2571_14.fastq
│   ├── FAP70729_pass_barcode11_003e2571_15.fastq
│   ├── FAP70729_pass_barcode11_003e2571_1.fastq
│   ├── FAP70729_pass_barcode11_003e2571_2.fastq
│   ├── FAP70729_pass_barcode11_003e2571_3.fastq
│   ├── FAP70729_pass_barcode11_003e2571_4.fastq
│   ├── FAP70729_pass_barcode11_003e2571_5.fastq
│   ├── FAP70729_pass_barcode11_003e2571_6.fastq
│   ├── FAP70729_pass_barcode11_003e2571_7.fastq
│   ├── FAP70729_pass_barcode11_003e2571_8.fastq
│   └── FAP70729_pass_barcode11_003e2571_9.fastq
└── reference
    ├── GCF_000001635.27_GRCm39_genomic.fna
    ├── GCF_000001635.27_GRCm39_genomic.fna.amb
    ├── GCF_000001635.27_GRCm39_genomic.fna.ann
    ├── GCF_000001635.27_GRCm39_genomic.fna.bwt
    ├── GCF_000001635.27_GRCm39_genomic.fna.fai
    ├── GCF_000001635.27_GRCm39_genomic.fna.pac
    └── GCF_000001635.27_GRCm39_genomic.fna.sa
```

### Outputs:

The output structure is as follows:

```
outputs/
├── config.yaml
├── mapped
│   ├── MAPPED_READ_ID_1
│   ├── ................
│   └── MAPPED_READ_ID_N
├── merged.bam
├── merged_reformatted.bam
├── merged_sorted.bam
├── stats
│   ├── coverage_table.txt
│   ├── custom_coverage_table.csv
│   ├── merged_sorted_stats.txt
│   └── samtools_coverage.txt
└── tables
    ├── alignment_table_chrom.csv
    ├── alignment_table.csv
    ├── digested_fragments_table.csv
    ├── filtered_alignment_table.csv
    ├── paohviz_output.csv
    └── raw_alignment_table.csv
```

