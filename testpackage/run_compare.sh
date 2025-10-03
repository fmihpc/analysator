#!/bin/bash -l



folder_1="/wrk-vakka/turso/group/spacephysics/CI_analysator/analysator_testpackage/verification_sets/380989e7fa7a331fb90c7ac6c496ebec6397dec9/"
folder_2="${PWD}/produced_plots/"
echo $(ls $(folder_2))
output=$((python3 ../analysator/testpackage/testpackage_compare.py ${folder_1} ${folder_2}) 2>&1)


if [[ $output == "" ]]; then
    echo "No differences found"
    exit 0
else
    echo "$output"
    exit 1
fi

