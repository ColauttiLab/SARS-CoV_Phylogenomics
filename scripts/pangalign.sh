#!/bin/bash
#SBATCH --mem=234G
#SBATCH -t 24:00:00
#SBATCH -c 24


module load nixpkgs/16.09 gcc/7.3.0
module load nixpkgs/16.09 intel/2018.3
module load mafft/7.397

date

# run mafftalign.sh FIRSTto generate ./intermediatedata/WuhanIM.afa

mafft --thread -1 --anysymbol --keeplength --add ./intermediatedata/WuhanIM.afa ./inputdata/pangolin/anonymised.encrypted.aln.putative.fasta > ./intermediatedata/WuhanPIM.afa

date
