#!/bin/bash

# run assigner
rna=$1
./scripts/assigner.R -i 10 -f 5 -s 0.75 -o data/test_assigned data/larmord_${rna}.txt data/observed_shifts_${rna}.txt data/larmord_accuracy_resname_nucleus.txt 