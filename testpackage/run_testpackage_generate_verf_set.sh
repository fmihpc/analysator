#!/bin/bash -l
#SBATCH -t 00:60:00
#SBATCH -J analysator_testpackage
#SBATCH -p short
#SBATCH -n 1
#SBATCH --array=0-10
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=16000

jobcount=$(( $SLURM_ARRAY_TASK_MAX - $SLURM_ARRAY_TASK_MIN + 1 )) 
index=$(( $SLURM_ARRAY_TASK_ID - $SLURM_ARRAY_TASK_MIN ))

hostname

source CI_env/bin/activate

export PTNONINTERACTIVE=1
export PTOUTPUTDIR=$1

python testpackage_colormap.py $jobcount $index
python testpackage_vdf.py $jobcount $index

echo Job $SLURM_ARRAY_TASK_ID complete.
echo "EXIT_CODE_FROM_JOB $?"
