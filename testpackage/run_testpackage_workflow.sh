#!/bin/bash -l
#SBATCH -t 00:45:00
#SBATCH -J analysator_testpackage
#SBATCH --constraint="carrington|ukko"
#SBATCH -p short
#SBATCH -n 1
#SBATCH --array=1-10
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=64000


#THIS SHOULD ONLY BE USED FOR GITHUB WORKFLOW TESTS
jobcount=$(( $SLURM_ARRAY_TASK_MAX - $SLURM_ARRAY_TASK_MIN + 1 )) 
index=$(( $SLURM_ARRAY_TASK_ID - $SLURM_ARRAY_TASK_MIN ))


hostname


module purge            
export PATH=/wrk-vakka/group/spacephysics/proj/appl/tex-basic/texlive/2023/bin/x86_64-linux:$PATH
older_python=false
echo "SLURM_JOB_ID=$SLURM_ARRAY_JOB_ID" >> $GITHUB_OUTPUT
if [[ $1 == "old_python" ]]; then
  module load $2
  older_python=true
else
  module load Python/3.10.4-GCCcore-11.3.0
fi
source CI_env/bin/activate

module list

mkdir -p $PWD/produced_plots/

export PTNONINTERACTIVE=1
export PTOUTPUTDIR=$PWD/produced_plots/


if $older_python; then
  echo "::warning:: Running with older python version $2"
  python ./testpackage/testpackage_commons.py $jobcount $index 
  echo "EXIT_CODE_FROM_JOB $?"
  exit 0
fi
python ./testpackage/testpackage_commons.py $jobcount $index $@
echo "EXIT_CODE_FROM_JOB $?"
