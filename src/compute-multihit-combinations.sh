#!/bin/bash
start=$(date +%s)

cancer=$1    
echo "Collecting data for " $cancer
mkdir ../result/$cancer
./acc2hit ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.training ../data/maf2dat-moderate/manifest_normal_normal.txt.training.txt.geneSampleList .1 > ../result/$cancer/$cancer-combinations
## Uses the geneSampleList input (not the BAM manifest) to align normals correctly!
./validate2hit ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.test ../result/$cancer/$cancer-combinations
end=$(date +%s)
echo "Total runtime: $((end - start)) seconds"
