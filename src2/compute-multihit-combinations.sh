#!/bin/bash
start=$(date +%s)

cancer=$1
serial=$2
echo "Collecting data for " $cancer
mkdir ../result/$cancer
if [[ "$serial" == "serial" ]]; then
echo "Doing in serial" $cancer
./acc2hitSerial ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.training ../data/maf2dat-moderate/manifest_normal_normal.txt.training.txt.geneSampleList .1 > ../result/$cancer/$cancer-combinations-serial
./validate2hit ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.test ../result/$cancer/$cancer-combinations-serial
python3 validateNormal2hit.py ../data/maf2dat-moderate/manifest_normal_normal.txt.test.txt.geneSampleList ../result/$cancer/$cancer-combinations-serial
elif [ "$type" == "sparse" ]; then
echo "Doing with sparse matrix - SIMD" $cancer
./sparseacc2hit ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.training ../data/maf2dat-moderate/manifest_normal_normal.txt.training.txt.geneSampleList .1 > ../result/$cancer/$cancer-combinations-sparse-simd
./validate2hit ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.test ../result/$cancer/$cancer-combinations-sparse-simd
python3 validateNormal2hit.py ../data/maf2dat-moderate/manifest_normal_normal.txt.test.txt.geneSampleList ../result/$cancer/$cancer-combinations-sparse-simd
else
echo "Doing with dense matrix - SIMD" $cancer
./acc2hit ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.training ../data/maf2dat-moderate/manifest_normal_normal.txt.training.txt.geneSampleList .1 > ../result/$cancer/$cancer-combinations-simd
./validate2hit ../data/maf2dat-moderate/$cancer.maf2dat.matrix.out.test ../result/$cancer/$cancer-combinations-dense-simd
python3 validateNormal2hit.py ../data/maf2dat-moderate/manifest_normal_normal.txt.test.txt.geneSampleList ../result/$cancer/$cancer-combinations-dense-simd
fi
# Uses the geneSampleList input (not the BAM manifest) to align normals correctly!

end=$(date +%s)
echo "Total runtime: $((end - start)) seconds"
