#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --job-name=FHA
#SBATCH --partition=short
#SBATCH -M carrington
##SBATCH -A spacephysics
#SBATCH --nodes=1
#SBATCH -c 1                  # CPU cores per task
#SBATCH -n 1                  # number of tasks (4xnodes)
#SBATCH --mem=12G # memory per node
#SBATCH --array=0-199
#SBATCH --mail-type=ALL

frameStart=1001  # set to initial frame
frameEnd=1612 # set to the final frame
#frameStart=1300  # set to initial frame
#frameEnd=1598 # set to the final frame
root=bulk1

# How many jobs? SLURM_ARRAY_TASK_COUNT does not work on all systems
# so calculate job count (or set it manually to match the array
# argument given to sbatch).
jobcount=$(( $SLURM_ARRAY_TASK_MAX - $SLURM_ARRAY_TASK_MIN + 1 )) 

# find job array index
index=$(( $SLURM_ARRAY_TASK_ID - $SLURM_ARRAY_TASK_MIN ))

frameEndC=$(( $frameEnd + 1 )) # Need to iterate to 1 past final frame
totalFrames=$(( $frameEndC - $frameStart )) # Total frames to calculate
increment=$(( $totalFrames / $jobcount )) # amount of frames per job (rounded down)

# Calculate remainder
remainder=$(( $totalFrames - $jobcount * $increment ))

start=$(( $frameStart + $index * $increment ))
end=$(( $start + $increment ))

# Remainder frames are divvied out evenly among tasks
if [ $index -lt $remainder ];
then 
    start=$(( $start + $index ))
    end=$(( $end + $index + 1 ))
else
    start=$(( $start + $remainder ))
    end=$(( $end + $remainder ))
fi;

# Ensure final job gets correct last frame
if [ $SLURM_ARRAY_TASK_ID -eq $SLURM_ARRAY_TASK_MAX ];
then 
    echo Verifying final frame: $end $frameEndC
    end=$frameEndC
fi;


module purge                 ## Purge modules for a clean start
module load Python/3.7.2-GCCcore-8.2.0
module load matplotlib
export PTNOLATEX=1 PTBACKEND=Agg PTNONINTERACTIVE=1 PYTHONPATH=$PYTHONPATH:/proj/ykempf/analysator/:/proj/ykempf/.local PYTHONUSERBASE=/proj/ykempf/.local
#export PTNOLATEX=1 PTBACKEND=Agg PTNONINTERACTIVE=1 PYTHONPATH=$PYTHONPATH:/proj/ykempf/analysator_master/:/proj/ykempf/.local PYTHONUSERBASE=/proj/ykempf/.local
umask 007


python plot_multi_with_xo.py $root $start $end 7
