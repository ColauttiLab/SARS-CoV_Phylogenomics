#!/bin/bash
#SBATCH --mem=60G
#SBATCH -t 24:00:00
#SBATCH -c 16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rc91@queensu.ca

# Memory usage ~ L^2 (square of sequence length)
# 30kb^2 ~ 8Gb

module load nixpkgs/16.09
module load gcc/7.3.0
module load blast+/2.10.1

date

sed 's/>2019-nCoV_MN[0-9]\+.*Sample \([0-9]\+\).*/>S\1i/g' < ./inputdata/IonTorrent_MinION.fasta > ./inputdata/IonTorrent_MinION.fasta.fix
sed -i 's/Sample_\([0-9]\+\).*/>S\1m/g' ./inputdata/IonTorrent_MinION.fasta.fix
sed -i 's/-/N/g' ./inputdata/IonTorrent_MinION.fasta.fix

sed 's/>.*EPI_ISL_//g' < ./inputdata/msa_0901/msa_0901.fasta > ./inputdata/msa_0901/msa_0901.fasta.fix
sed -i 's/\([0-9]\+\).*/>\1/g' ./inputdata/msa_0901/msa_0901.fasta.fix
sed -i 's/-/N/g' ./inputdata/msa_0901/msa_0901.fasta.fix

date

blastn -query ./inputdata/IonTorrent_MinION.fasta.fix\
 -subject ./inputdata/msa_0901/msa_0901.fasta.fix\
 -strand plus -word_size 1000 -gapopen 0 -gapextend 0\
 -max_hsps 10000 -max_target_seqs 10000\
 -out BlastResults.out -outfmt 6

date
