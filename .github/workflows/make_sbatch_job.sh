#!/usr/bin/env bash
if [[ ! $ARRAY_SIZE ]]; then
    ARRAY_SIZE=14 
fi
sbatch -W -o --array=1-$ARRAY_SIZE "$1" ./testpackage/run_compare.sh old_python ${{ steps.pyversion.outputs.PYTHON }} > jobid_$1 || srun --pty cat $1 || cat $1
srun --pty cat $1 || cat $1
export JOBID=$(srun grep -Po '\d+' jobid1.txt || grep -Po '\d+' jobid_$1)
export SACCT_LOG=$(sacct -j $JOBID -o job,state,node | grep FAILED) 
if [[ $SACCT_LOG ]]; then
  echo "Some job failed on a node, try to take a look at the slurm log."
  echo "$SACCT_LOG" 
  exit 1
fi
