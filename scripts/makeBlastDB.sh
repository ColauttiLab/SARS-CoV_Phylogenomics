#!/bin/bash
makeblastdb -in ./inputdata/nextstrain_sequences.fasta -parse_seqids -dbtype 'nucl' -hash_index -title 'nextstrainDB' -out './intermediatedata/nextstrain_sequences/nextstrainDB'
