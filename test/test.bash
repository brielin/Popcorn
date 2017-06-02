#!/bin/bash
# basic compute commands
python ../__main__.py compute --bfile EAS_ref EAS.cscore
python ../__main__.py compute --bfile1 EAS_ref --bfile2 EUR_ref EAS_EUR.cscore

# basic fit commands
python ../__main__.py fit --sfile EAS_test.txt --cfile EAS.cscore EAS.out
python ../__main__.py fit --sfile1 EAS_test.txt --sfile2 EUR_test.txt --cfile EAS_EUR.cscore EAS_EUR.out

# exmples illustrating some options
python ../__main__.py compute --bfile1 EAS_ref --bfile2 EUR_ref --maf 0.05 --gen_effect EAS_EUR.cscore
python ../__main__.py fit --sfile1 EAS_test.txt --sfile2 EUR_test.txt --maf 0.05 --gen_effect \
    --K1 0.1 --K2 0.2 --P1 0.5 --P2 0.4 --cfile EAS_EUR.cscore EAS_EUR.out
