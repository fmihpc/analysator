#!/bin/bash -l



#verification set is from: 6.10.2025 (vdf 15.10.2025 from image_compare branch)
folder_1="/wrk-vakka/turso/group/spacephysics/CI_analysator/analysator_testpackage/verification_sets/380989e7fa7a331fb90c7ac6c496ebec6397dec9/$1/"
folder_2="${PWD}/produced_plots/$1/"
python3 ../analysator/testpackage/testpackage_compare.py ${folder_1} ${folder_2} && echo "No differences found"
