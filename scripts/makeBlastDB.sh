#!/bin/bash
makeblastdb -in ./inputdata/nextstrain_sequences.fasta -parse_seqids -dbtype 'nucl' -hash_index -title "nextstraindb" 
