#!/bin/bash -l
#SBATCH -t 00:30:00
#SBATCH -J analysator_testpackage
#SBATCH -p short
#SBATCH -n 1
#SBATCH --array=1-2
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=16000


#THIS SHOULD ONLY BE USED FOR GITHUB WORKFLOW TESTS
jobcount=$(( $SLURM_ARRAY_TASK_MAX - $SLURM_ARRAY_TASK_MIN + 1 )) 
index=$(( $SLURM_ARRAY_TASK_ID - $SLURM_ARRAY_TASK_MIN ))


hostname


module purge            
export PATH=/wrk-vakka/group/spacephysics/proj/appl/tex-basic/texlive/2023/bin/x86_64-linux:$PATH
module load Python/3.10.4-GCCcore-11.3.0
module load ImageMagick/7.1.0-37-GCCcore-11.3.0



source CI_env/bin/activate

module list

mkdir -p $PWD/produced_plots

export PTNONINTERACTIVE=1
export PTOUTPUTDIR=$PWD/produced_plots/


python ./testpackage/color_small_test.py $jobcount $index
echo "EXIT_CODE_FROM_JOB $?"
#echo Job $SLURM_ARRAY_TASK_ID complete.
