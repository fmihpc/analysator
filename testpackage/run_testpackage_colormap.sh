#!/bin/bash -l
#SBATCH -t 00:20:00
#SBATCH -J analysator_testpackage
#SBATCH -p serial
#SBATCH -n 1
#SBATCH --array=0-99
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=16000
#SBATCH --constraint=hsw

jobcount=$(( $SLURM_ARRAY_TASK_MAX - $SLURM_ARRAY_TASK_MIN + 1 )) 
index=$(( $SLURM_ARRAY_TASK_ID - $SLURM_ARRAY_TASK_MIN ))

#module load mayavi2
module purge
module load python-env/3.5.3

export PTNONINTERACTIVE=1
python testpackage_colormap.py $jobcount $index
echo Job $SLURM_ARRAY_TASK_ID complete.
