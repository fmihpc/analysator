#!/bin/bash -l



#verification set is from: 6.10.2025 
folder_1="/wrk-vakka/turso/group/spacephysics/CI_analysator/analysator_testpackage/verification_sets/e39b15ed638dee8235f87d88cb12221f89d2e9d7/"
folder_2="${PWD}/produced_plots/"
output=$((python3 ../analysator/testpackage/testpackage_compare.py ${folder_1} ${folder_2}) 2>&1)


if [[ $output == "" ]]; then
    echo "No differences found"
    exit 0
else
    echo "$output"
    exit 1
fi

