#!/bin/bash
#SBATCH --mem=234G
#SBATCH -t 72:00:00
#SBATCH -c 24
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rc91@queensu.ca

date

module load gcc/8.3.0
module load r/4.0.0

Rscript --vanilla Distcalc.R

date
