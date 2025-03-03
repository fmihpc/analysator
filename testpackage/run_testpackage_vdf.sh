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

#module load mayavi2
# module purge
# module load python-env/3.5.3

source ~/pyvenv

export PTNONINTERACTIVE=1
export PTOUTPUTDIR=/wrk-vakka/users/mjalho/analysator_testpackage/

python testpackage_vdf.py $jobcount $index
echo Job $SLURM_ARRAY_TASK_ID complete.
