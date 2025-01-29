#!/bin/bash
#SBATCH -c 64
#SBATCH --mem 720GB
#SBATCH -p long,big-mem,normal,express

source ~/.bashrc
conda activate notebook

snakemake --profile profile/