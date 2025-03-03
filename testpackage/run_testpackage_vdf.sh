#!/bin/bash -l
#SBATCH -t 00:20:00
#SBATCH -J analysator_testpackage_vdf
#SBATCH -p short
#SBATCH -n 1
#SBATCH --array=1
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=16000

jobcount=$(( $SLURM_ARRAY_TASK_MAX - $SLURM_ARRAY_TASK_MIN + 1 )) 
index=$(( $SLURM_ARRAY_TASK_ID - $SLURM_ARRAY_TASK_MIN ))

source /wrk-vakka/group/spacephysics/proj/analysator_testpackage/pyvenv.sh

export PTNONINTERACTIVE=1
export PTOUTPUTDIR=$PWD/

python testpackage_vdf.py $jobcount $index
echo Job $SLURM_ARRAY_TASK_ID complete.
