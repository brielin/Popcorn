#!/bin/bash
# basic compute commands
popcorn compute --bfile EAS_ref EAS.cscore
popcorn compute --bfile1 EAS_ref --bfile2 EUR_ref EAS_EUR.cscore

# basic fit commands
popcorn fit --sfile EAS_test.txt --cfile EAS.cscore EAS.out
popcorn fit --sfile1 EAS_test.txt --sfile2 EUR_test.txt --cfile EAS_EUR.cscore EAS_EUR.out

# exmples illustrating some options
popcorn compute --bfile1 EAS_ref --bfile2 EUR_ref --maf 0.05 --gen_effect EAS_EUR.cscore
popcorn fit --sfile1 EAS_test.txt --sfile2 EUR_test.txt --maf 0.05 --gen_effect \
    --K1 0.1 --K2 0.2 --P1 0.5 --P2 0.4 --cfile EAS_EUR.cscore EAS_EUR.out
