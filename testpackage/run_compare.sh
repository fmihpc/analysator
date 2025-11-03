#!/bin/bash -l


verf_loc="/wrk-vakka/turso/group/spacephysics/CI_analysator/analysator_testpackage/verification_sets"

#gets latest verfication set (based on modification date -> grep directories only -> take firstline -> get last word)
folder_1="$verf_loc/$(ls -lth $verf_loc | grep ^d | head -n1 | grep -Po '\w+$')/$1/" 
folder_2="${PWD}/produced_plots/$1/"
python3 ../analysator/testpackage/testpackage_compare.py ${folder_1} ${folder_2} && echo "No differences found in produced images"
