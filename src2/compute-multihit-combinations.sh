#!/bin/bash
start=$(date +%s)

cancer=$1
serial=$2
echo "Collecting data for " $cancer
mkdir ../result/$cancer
if [[ "$serial" == "serial" ]]; then
echo "Done in serial" $cancer
./acc2hitSerial ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.training ../data/maf2dat-moderate/manifest_normal_normal.txt.training.txt.geneSampleList .1 > ../result/$cancer/$cancer-combinations-serial
./validate2hit ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.test ../result/$cancer/$cancer-combinations-serial
python3 validateNormal2hit.py ../data/maf2dat-moderate/manifest_normal_normal.txt.test.txt.geneSampleList ../result/$cancer/$cancer-combinations-serial
else
echo "Done in parallel" $cancer
./acc2hit ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.training ../data/maf2dat-moderate/manifest_normal_normal.txt.training.txt.geneSampleList .1 > ../result/$cancer/$cancer-combinations-simd
./validate2hit ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.test ../result/$cancer/$cancer-combinations-simd
python3 validateNormal2hit.py ../data/maf2dat-moderate/manifest_normal_normal.txt.test.txt.geneSampleList ../result/$cancer/$cancer-combinations-simd
fi
# Uses the geneSampleList input (not the BAM manifest) to align normals correctly!

end=$(date +%s)
echo "Total runtime: $((end - start)) seconds"
