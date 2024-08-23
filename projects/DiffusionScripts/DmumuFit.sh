#!/bin/bash -l
#SBATCH --job-name=$jobName$
#SBATCH --time=00:10:00
#SBATCH --account=project_2000203
#SBATCH --partition=medium
#SBATCH --nodes=1 ## number of nodes

module load python-data

# set the number of threads based on --cpus-per-task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Bind OpenMP threads to hardware threads
export OMP_PLACES=cores

#echo Calculating a total of $totalFrames frames divided amongst $jobcount jobs.
#echo Using a remainder of $remainder.
#echo Current job id is $SLURM_ARRAY_TASK_ID and calculated frames are
#echo from $start to $end 

export PTNONINTERACTIVE=1

python3 /users/dubartma/analysis/subgrid/TimeDependence/solverIvascenko.py dmumu $folder$ $start$ $end$ $interval$ $twindow$

echo Job $SLURM_ARRAY_TASK_ID complete.
~                                       
