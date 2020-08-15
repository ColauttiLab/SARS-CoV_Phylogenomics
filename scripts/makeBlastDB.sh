#!/bin/bash
makeblastdb -in sequences_0815.fasta -parse_seqids -dbtype 'nucl' -hash_index -title 'nextstrainDB' -out 'nextstrainDB'
