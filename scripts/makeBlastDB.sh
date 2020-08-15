#!/bin/bash
makeblastdb -in nextstrain_sequences_051520.fasta -parse_seqids -dbtype 'nucl' -hash_index -title "nextstraindb05152020" 
