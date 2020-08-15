#!/bin/bash
FILE="outputs/seqs.txt"
if [ ! -f "$FILE" ]
then
    touch $FILE
elif [ -s "$FILE" ]
then
    > "$FILE"
fi

# NOTE blast directory may vary by user

while read -r CURRENT_LINE
    do
        /usr/local/ncbi/blast/bin/blastdbcmd -db ./inputdata/nextstrain_sequences_051520/nextstrain_sequences_051520.fasta -entry $CURRENT_LINE | tee -a ./seqs.txt
        ((LINE++))
done < "outputs/GBACC.txt"
