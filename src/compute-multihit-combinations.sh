#!/bin/bash
cancer=$1    
echo "Collecting data for " $cancer
mkdir ../result/$cancer
./acc2hit ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.training ../data/maf2dat-moderate/manifest_normal_normal.txt.training.txt.geneSampleList .1 > ../result/$cancer/$cancer-combinations
./validate2hit.sh ../result/$cancer/$cancer-combinations $cancer

