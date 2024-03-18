#!/bin/bash
#SBATCH --job-name=$jobName$
#SBATCH --time=15:00:00
#SBATCH --account=project_2000203
#SBATCH --partition=medium
#SBATCH --nodes=8             ## Number of nodes
#SBATCH --ntasks-per-node=16   ## number of tasks 
#SBATCH --cpus-per-task=16    ## CPU cores per task
#SBATCH --hint=multithread

# set the number of threads based on --cpus-per-task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Bind OpenMP threads to hardware threads
export OMP_PLACES=cores

module load boost papi jemalloc
module load python-data

umask 007
cd ./bulk/

executable="/users/dubartma/vlasiator/vlasiator"
configfile=$configfile$

srun $executable --run_config $configfile

#echo "setting a dependent job for task $SLURM_JOB_ID"
#sbatch --dependency=afterany:$SLURM_JOB_ID $variables$


