#!/bin/bash -l
#SBATCH -t 00:60:00
#SBATCH -J analysator_testpackage
#SBATCH -p short
#SBATCH -n 1
#SBATCH --array=0-10
#SBATCH --constraint="carrington|ukko"
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=16000

jobcount=$(( $SLURM_ARRAY_TASK_MAX - $SLURM_ARRAY_TASK_MIN + 1 )) 
index=$(( $SLURM_ARRAY_TASK_ID - $SLURM_ARRAY_TASK_MIN ))

hostname

echo "SLURM_JOB_ID=$SLURM_ARRAY_JOB_ID" >> $GITHUB_OUTPUT
source CI_env/bin/activate
export PATH=/wrk-vakka/group/spacephysics/proj/appl/tex-basic/texlive/2023/bin/x86_64-linux:$PATH

export PTNONINTERACTIVE=1
export PTOUTPUTDIR=$1

python ./testpackage/testpackage_commons.py $jobcount $index

exit $?

