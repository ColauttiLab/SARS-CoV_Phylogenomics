#!/bin/bash
sed 's/-/N/g' ./inputdata/msa_0901/msa_0901.fasta > ./intermediatedata/msa_Ns.fasta
makeblastdb -in ./intermediatedata/msa_Ns.fasta -dbtype 'nucl' -hash_index -title 'nextstrainDB' -out './intermediatedata/nextstrain_sequences/nextstrainDB'
