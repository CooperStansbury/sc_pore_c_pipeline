#!/bin/bash

##################
####  Slurm preamble

#### #### ####  These are the most frequently changing options

####  Request resources here
####    These are typically, number of processors, amount of memory,
####    an the amount of time a job requires.  May include processor
####    type, too.

#SBATCH --nodes=36
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=180GB
#SBATCH --time=24:00:00

####  Slurm account and partition specification here
####    These will change if you work on multiple projects, or need
####    special hardware, like large memory nodes or GPUs.

#SBATCH --account=indikar1
#SBATCH --partition=standard

#### #### ####  These are the least frequently changing options

####  Your e-mail address and when you want e-mail

#SBATCH --mail-user=cstansbu@umich.edu
#SBATCH --mail-type=END,FAIL

#### where to write log files

#SBATCH --output=log/hpc/slurm-%j_%x.out
snakemake --profile config/slurm --latency-wait 90 