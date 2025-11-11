#!/bin/bash -l

verf_loc="/wrk-vakka/turso/group/spacephysics/CI_analysator/analysator_testpackage/verification_sets"

#if pass we do not check for anything
if echo $@ | grep -q -P "\spass$|\spass\s"; then
   exit 0
fi

check=true

for i in $@
do
    check=false
    echo "Comparing for $i"
    #gets latest verfication set (based on modification date -> grep directories only -> take firstline -> get last word)
    folder_1="$verf_loc/$(ls -lth $verf_loc | grep ^d | head -n1 | grep -Po '\w+$')/testpackage_run/$i/" 
    folder_2="${PWD}/produced_plots/testpackage_run/$i/"
    python3 ../analysator/testpackage/testpackage_compare.py ${folder_1} ${folder_2} && echo "No differences found in produced images"
done

if $check;
then
    echo "Comparing all"
    folder_1="$verf_loc/$(ls -lth $verf_loc | grep ^d | head -n1 | grep -Po '\w+$')/testpackage_run/" 
    folder_2="${PWD}/produced_plots/testpackage_run/"
    python3 ../analysator/testpackage/testpackage_compare.py ${folder_1} ${folder_2} && echo "No differences found in produced images"
fi
