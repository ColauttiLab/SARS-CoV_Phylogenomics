#!/bin/bash
FILE="./outputs/seqs.txt"
if [ ! -f "$FILE" ]
then
    touch $FILE
elif [ -s "$FILE" ]
then
    > "$FILE"
fi

while read -r CURRENT_LINE
    do
        /usr/local/ncbi/blast/bin/blastdbcmd -db ./intermediatedata/nextstrain_sequences/nextstrainDB -entry $CURRENT_LINE | tee -a ./outputs/seqs.txt
        ((LINE++))
done < "./outputs/GBACC.txt"
