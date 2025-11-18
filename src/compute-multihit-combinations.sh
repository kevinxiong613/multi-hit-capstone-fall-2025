#!/bin/bash
cancer=$1    
echo "Collecting data for " $cancer
mkdir ../result/$cancer
# Calls it with training data, normal data, some beta value 0.1, and outputs to a file
./acc2hit ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.training ../data/maf2dat-moderate/manifest_normal_normal.txt.training.txt.geneSampleList .1 > ../result/$cancer/$cancer-combinations

# I think this is wrong? In code it says it needs test data and combinations file
# ./validate2hit ../result/$cancer/$cancer-combinations $cancer
./validate2hit ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.test ../result/$cancer/$cancer-combinations