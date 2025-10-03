#!/bin/bash -l



folder_1="/wrk-vakka/turso/group/spacephysics/CI_analysator/analysator_testpackage/verification_sets/380989e7fa7a331fb90c7ac6c496ebec6397dec9/"
folder_2="/wrk-vakka/turso/group/spacephysics/CI_analysator/analysator_testpackage/verification_sets/5e40f1f7621c858984903d89771a5e7a00047863/"


output=$(python3 ../analysator/testpackage/testpackage_compare.py ${folder_1} ${folder_2} )
#Make it so stderr is also captured


if [[ $output == "" ]]; then
    echo "No differences found"
    exit 0
else
    echo "$output"
    exit 1
fi

