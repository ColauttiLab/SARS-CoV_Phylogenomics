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
        /usr/local/ncbi/blast/bin/blastdbcmd -db ~/ColauttiLabScratch/COVID-19/August/QGLO_COVID/inputdata/nextstrainDB/nextstrainDB -entry $CURRENT_LINE | tee -a outputs/seqs.txt
        ((LINE++))
done < "outputs/GBACC.txt"
