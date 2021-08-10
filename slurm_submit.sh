#!/bin/bash

#SBATCH --mail-user=cstansbu@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180gb
#SBATCH --time=00-24:00:00
#SBATCH --account=indikar1
#SBATCH --partition=standard
#SBATCH --output=/home/cstansbu/git_repositories/sc_pore_c_pipeline/logs/%x-%j.logs

snakemake --profile config/slurm