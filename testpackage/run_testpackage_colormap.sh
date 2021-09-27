#!/bin/bash -l
#SBATCH -t 00:30:00
#SBATCH -J analysator_testpackage
#SBATCH -p short
#SBATCH -n 1
#SBATCH --array=0-40
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=16000

jobcount=$(( $SLURM_ARRAY_TASK_MAX - $SLURM_ARRAY_TASK_MIN + 1 )) 
index=$(( $SLURM_ARRAY_TASK_ID - $SLURM_ARRAY_TASK_MIN ))

#module load mayavi2
module purge
#module load python-env/3.5.3
module load Python/3.7.2-GCCcore-8.2.0

export PATH=/proj/jesuni/projappl/tex-basic/texlive/2020/bin/x86_64-linux:$PATH
module load matplotlib

export PTNONINTERACTIVE=1
#export PTNOLATEX=1
export PTNOLATEX=
export PTOUTPUTDIR=/wrk/users/markusb/Plots/

python testpackage_colormap.py $jobcount $index
echo Job $SLURM_ARRAY_TASK_ID complete.
