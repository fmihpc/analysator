#!/bin/bash -l
#SBATCH -t 00:30:00
#SBATCH -J analysator_testpackage_compare
#SBATCH --constraint="ukko|carrington"
#SBATCH -p short
#SBATCH -n 1
#SBATCH --array=1-10
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=16000


#THIS SHOULD ONLY BE USED FOR GITHUB WORKFLOW TESTS
jobcount=$(( $SLURM_ARRAY_TASK_MAX - $SLURM_ARRAY_TASK_MIN + 1 )) 
index=$(( $SLURM_ARRAY_TASK_ID - $SLURM_ARRAY_TASK_MIN ))


module purge
module load Python/3.10.4-GCCcore-11.3.0
source CI_env/bin/activate
module load libglvnd/1.7.0-GCCcore-13.3.0
module list

verf_loc="/wrk-vakka/turso/group/spacephysics/CI_analysator/analysator_testpackage/verification_sets"

#if pass we do not check for anything
if echo $@ | grep -q -P "\spass$|\spass\s"; then
   exit 0
fi

check=true

#Note that this is skipped if no arguments are passed
for i in $@
do
    check=false
    echo "Comparing for $i"
    #gets latest verfication set (based on modification date -> grep directories only -> take firstline -> get last word)
    folder_1="$verf_loc/$(ls -lth $verf_loc | grep ^d | head -n1 | grep -Po '\w+$')/$i/" 
    folder_2="${PWD}/produced_plots/$i/"
    python3 ../analysator/testpackage/testpackage_compare.py ${folder_1} ${folder_2} $jobcount $index && echo "No differences found in produced images"
done

if $check;
then
    echo "Comparing all"
    folder_1="$verf_loc/$(ls -lth $verf_loc | grep ^d | head -n1 | grep -Po '\w+$')/" 
    folder_2="${PWD}/produced_plots/"
    python3 ../analysator/testpackage/testpackage_compare.py ${folder_1} ${folder_2} $jobcount $index && echo "No differences found in produced images"
fi
