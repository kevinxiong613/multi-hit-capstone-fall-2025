## Overview
* Cancer is known to result from a combination of a small number of genetic defects. 
* The specific combinations of mutations responsible for the vast majority of cancers have not been identified. 
* We present a fundamentally different approach for identifying the cause of individual instances of cancer: we search for combinations of carcinogenic mutations (multi-hit combinations) instead of driver mutations. 
* We developed an algorithm that identified a set of multi-hit combinations that differentiate between tumor and normal tissue samples with $91\%$ sensitivity (95% Confidence Interval (CI)=89-92%) and 93%  specificity (95% CI=91-94%) on average for seventeen cancer types. 

Note: If doing on rlogin or glogin, change src/Makefile to be "gcc" instead of "CC = gcc-15". CC = gcc-15 was only used to be compatible with our Macbooks.

## Identify Combinations for One Cancer Type - Serial performance

If on Mac OS, do:
```
brew install gcc
```
Otherwise if you're in a linux environment with gcc, change the makefile to have 
```
gcc
```
instead of
```
gcc-15
```

```
cd src 
make
./compute-multihit-combinations.sh BRCA serial
```

## Identify Combinations for One Cancer Type - Parallel performance with sparse matrix
```
cd src 
make
./compute-multihit-combinations.sh BRCA sparse
```

## Identify Combinations for One Cancer Type - Parallel performance with dense matrix
```
cd src 
make
./compute-multihit-combinations.sh BRCA [anything that isn't serial or sparse]
```

## Identify Combinations for One Cancer Type - OpenACC SIMD Parallelization

Note: SIMD parallelization was done on glogin. If you want to reproduce our results, please do so in glogin.

```
cd src2
make
./compute-multihit-combinations.sh BRCA [anything that isn't serial or sparse]
```

## Verify correctness using Bernie's script

```
cd src OR cd src2
python3 addParanthesis.py [combinations file]

python3 verifyAccuracy.py ../data/maf2dat-moderate/BRCA.maf2dat.matrix.out.test ../data/maf2dat-moderate/manifest_normal_normal.txt.test.txt.geneSampleList [one of the 3 combinations files from above]
```

Note that running add paranthesis will make the original validate script not work anymore.

Running addParanthesis on the same combinations file will cause it to no longer work with verifyAccuracy. Make sure to only run it once per a combinations file!

The goal is to maximize tumor sample coverage and minimize normal sample coverage.

---

## Output
The output combinations are in the folder result/BRCA

---

