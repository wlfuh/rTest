#!/bin/bash

# run assigner
#./assigner.R -i 5000 -s 0.75 data/larmord_2KOC.txt data/observed_shifts_2KOC.txt data/larmord_accuracy_resname_nucleus.txt -o data/out/test.txt
for out in $(ls ../data/observed_shifts_*.txt)
do
    id=$(echo $out| cut  -d'_' -f 3 | cut -d'.' -f 1)
    obs="data/$out"
    pred="data/larmord_$id.txt"
    acc="data/larmord_accuracy_resname_nucleus.txt"
    res="data/out/larmord_$id_out.txt"
    ./assigner.R -i 10 -s 0.75 $pred $obs $acc -o $res
done
