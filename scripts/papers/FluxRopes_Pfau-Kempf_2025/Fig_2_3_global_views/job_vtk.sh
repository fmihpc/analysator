#!/bin/bash -l
#SBATCH -t 1:0:0
#SBATCH -J array
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-1
#SBATCH --no-requeue
#SBATCH --mem=0
#SBATCH -M carrington
#SBATCH --nice=1042


frameStart=1600  # set to initial frame
frameEnd=1601 # set to the final frame

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

source prologue.sh

export OMP_NUM_THREADS=64

srun --mpi=pmix_v3 python extract_vtk.py -t $start -tracing 0
