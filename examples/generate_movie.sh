#!/bin/bash -l
#SBATCH -t 00:30:00
#SBATCH -J array_movie
#SBATCH -p serial
#SBATCH -n 1
#SBATCH --array=0-200
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=16000
#SBATCH --constraint=hsw

# This script can be used on taito to generate an array job, which renders multiple frames
# in order to e.g. make a movie.

frameStart=0  # set to initial frame
frameEnd=2708 # set to the final frame
jobNumber=200  # equal to the higher value in the array argument to sbatch

frameEndC=$(( $frameEnd + 1 )) # Need to iterate to 1 past final frame
jobsC=$(( $jobNumber + 1 )) # Actually this many jobs
totalFrames=$(( $frameEndC - $frameStart )) # Total frames to calculate
increment=$(( $totalFrames / $jobNumber )) # amount of frames per job (rounded down)

# Calculate remainder
remainder=$(( $totalFrames - $jobsC * $increment ))

start=$(( $frameStart + $SLURM_ARRAY_TASK_ID * $increment ))
end=$(( $start + $increment ))

# Remainder frames are divvied out evenly among tasks
if [ $SLURM_ARRAY_TASK_ID -lt $remainder ];
then 
    start=$(( $start + $SLURM_ARRAY_TASK_ID ))
    end=$(( $end + $SLURM_ARRAY_TASK_ID + 1 ))
else
    start=$(( $start + $remainder ))
    end=$(( $end + $remainder ))
fi;


# Ensure final job gets correct last frame
if [ $SLURM_ARRAY_TASK_ID -eq $jobNumber ];
then 
    echo Verifying final frame: $end $frameEndC
    end=$frameEndC
fi;

#echo Calculating a total of $totalFrames frames divided amongst $jobsC jobs.
#echo Using a remainder of $remainder.
#echo Current job id is $SLURM_ARRAY_TASK_ID and calculated frames are
#echo from $start to $end 

module load mayavi2
python generate_panel.py $start $end
echo Job $SLURM_ARRAY_TASK_ID complete.
