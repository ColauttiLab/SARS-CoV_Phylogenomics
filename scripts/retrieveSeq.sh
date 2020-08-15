#!/bin/bash
FILE="outputs/seqs.txt"
if [ ! -f "$FILE" ]
then
    touch $FILE
elif [ -s "$FILE" ]
then
    > "$FILE"
fi

while read -r CURRENT_LINE
    do
        /usr/local/ncbi/blast/bin/blastdbcmd -db ~/Dropbox/_Synced/PROJECTS/COVID_Sequencing/QGLO_COVID/inputdata/nextstrain_sequences_051520/nextstrain_sequences_051520.fasta -entry $CURRENT_LINE | tee -a ./seqs.txt
        ((LINE++))
done < "outputs/GBACC.txt"
