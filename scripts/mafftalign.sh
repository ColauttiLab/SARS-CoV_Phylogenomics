#!/bin/bash
#SBATCH --mem=234G
#SBATCH -t 24:00:00
#SBATCH -c 24

module load gcc/7.3.0
module load intel/2018.3
module load mafft/7.397

date

# Align to nextstrain alignment
mafft --thread -1 --keeplength --add ./inputdata/IonTorrent_consensus.fasta ./inputdata/Wuhan.afa > ./intermediatedata/WuhanI.afa
mafft --thread -1 --keeplength --add ./inputdata/minion_consensus.fasta ./intermediatedata/WuhanI.afa > ./intermediatedata/WuhanIM.afa

# Align to Reference sequences
mafft --thread -1 --keeplength --add ./intermediatedata/WuhanIM.afa ./inputdata/RefSeqs.fasta > ./intermediatedata/RefWIM.afa


echo "Header line numbers:"

grep -n '>' ./intermediatedata/WuhanIM.afa | head -3

echo "Use above to check line of reference"

echo "Cutting at line 518 using sed '1,518d'" # Cuts Wuhan sequence to avoid duplicates in final file

sed '1,518d' ./intermediatedata/WuhanIM.afa > ./intermediatedata/Samples.afa

cp ./inputdata/msa_0901/msa_0901.fasta ./intermediatedata/BRaligned.afa

cat ./intermediatedata/Samples.afa >> ./intermediatedata/BRaligned.afa # Append aligned to NextStrain

echo "Wuhan cut and Samples appended to ./intermediatedata/BRaligned.afa"

echo "CHECK to make sure 518 is the first line of Sample IDs minus 1"

head -1 ./intermediatedata/Samples.afa # Check to make sure Wuhan sequence was properly cut

echo "Check above is Sample ID (head -1 ./intermediate/Samples.afa)"

date
