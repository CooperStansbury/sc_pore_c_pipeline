# Single Cell Pore-C Pipeline

## Overview 

This tool is intended to create analysis files for single-cell pore-C reads.

## Installation

The pipeline's environment can be built by calling the following `conda` functions from inside the `sc_pore_c_pipeline` directory:

```
conda env create
conda activate pipeline
```

## Usage

This is a [snakemake tool](https://snakemake.readthedocs.io/en/stable/index.html). It is designed to take reads after basecalling through to several analysis files. 

### Running the Pipeline

The whole pipeline may be run from inside the `sc_pore_c_pipeline` directory like so:

```
snakemake -j <n_cores>
```

Individual steps may be run, as long as they do not contain wildcards in their input/output arguments:

```
snakemake -j <n_cores> <rule_name>
```

### Configuration

The file `sc_pore_c_pipeline/config/config.yaml` contains the following parameters that must be specified for each run of the pipeline:

- **reads:** An input directory containing _only_ files with the extension `.fastq` or `.gz` (for zipped .`fastq`) files. The file type is specified by the 'fastq_filetype' parameter listed below. The directory path MUST end with a "/", for example: `inputs/` but not `inputs`.
- **reference:** A path to an indexed, unzipped reference genome (`.fna`) file for the alignment tool of choice. Note that the reference must be indexed by the tool used to align the genomes in the parameter 'aligner' below.
- **assembly:** A path to an "assembly" file from the [The Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc)  for example: [https://www.ncbi.nlm.nih.gov/grc/human/data](https://www.ncbi.nlm.nih.gov/grc/human/data). Note that this file is formated as a `.csv` file with the column headers: `Chromosome,Total length (bp),GenBank accession,RefSeq accession`. There are a few examples in `sc_pore_c_pipeline/assembly_files/`. 
- **output:** A path to a directory for outputs. The pipeline will create the folder if it does not exist. Note that the pipeline will not regenerate files that already do exist, if their inputs have not changed between pipeline executions.
- **aligner:** A string formated function for alignment. Note that the refenerence genome must be unzipped and indexed by the tool used to align. Currently, the environment supports both `minimap2` and `bwa` alignment tools. Common examples of the function include:

```
(1) bwa mem
(2) bwa bwasw -b 5 -q 2 -r 1 -T 15 -z 10
(3) minimap2 -ax map-ont
```

- **threads:** An integer number for the number of threads. This will depend on the computational resources available at runtime. 
- **chromosome:** An integer number for subprocessing of a single chromosome. Can also be "all". This affects ONLY outouts for [Paohvis](https://www.aviz.fr/Research/Paohvis).
- **assembly_field:** A string that must be either `GenBank` or `RefSeq` exactly. This set the field that the `alignment_table.csv` file is mapped to on the assemlby, and is critical for creating outputs for [Paohvis](https://www.aviz.fr/Research/Paohvis).
- **fastq_filetype:** Can be either `.fastq` or `.gz`. Must match the structure of the files in the 'reads' parameter. MUST include the period.

### Example Input/Output Structure

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
│   ├── FAP70729_pass_barcode11_003e2571_0.bam
│   ├── FAP70729_pass_barcode11_003e2571_10.bam
│   ├── FAP70729_pass_barcode11_003e2571_11.bam
│   ├── FAP70729_pass_barcode11_003e2571_12.bam
│   ├── FAP70729_pass_barcode11_003e2571_13.bam
│   ├── FAP70729_pass_barcode11_003e2571_14.bam
│   ├── FAP70729_pass_barcode11_003e2571_15.bam
│   ├── FAP70729_pass_barcode11_003e2571_1.bam
│   ├── FAP70729_pass_barcode11_003e2571_2.bam
│   ├── FAP70729_pass_barcode11_003e2571_3.bam
│   ├── FAP70729_pass_barcode11_003e2571_4.bam
│   ├── FAP70729_pass_barcode11_003e2571_5.bam
│   ├── FAP70729_pass_barcode11_003e2571_6.bam
│   ├── FAP70729_pass_barcode11_003e2571_7.bam
│   ├── FAP70729_pass_barcode11_003e2571_8.bam
│   └── FAP70729_pass_barcode11_003e2571_9.bam
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

